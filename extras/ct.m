function varargout = ct (varargin)
% dc3dm convergence tests. The tests in this file are used to create figures
% like Figs. 3 and 4 in the SRL paper. This file is probably not very easy to
% use. See ex.m for step-by-step instructions.
%
% Roughly, the procedure to run these tests is as follows:
%     s = ct('cts_run_final', 5); % 5 is for test number 5.
% Execute in a shell all the commands that are emitted. When these are done,
% run
%     s = ct('cts_a1_init', s);
%     figure; ct('cts_a1_p1', s);
%     figure; ct('cts_a1_show_spiky', s);
% That was the procedure for one test. To run each of five tests for each of
% 9 component combinations, run
%     ss = ct('cts_all_setup');
% Then do the above for each i: s = ss{i}.
  [varargout{1:nargout}] = feval(varargin{:});
end

% ------------------------------------------------------------------------------
% Hopefully only this section's code needs to be modified to run on your system.

function addpaths ()
  addpath ../matlab;
end

function e = get_env ()
  e.save_dir = '/scratch/ambrad/dctf1/cts';
  e.bin_dc3dm = '../bin/dc3dm';
end

% ------------------------------------------------------------------------------
% Convergence test of EIGA and AIGA for a smooth problem with a small mesh
% height, and various BCs.

function ss = cts_all_setup ()
% Run all the tests on all source-receiver combinations.
  e = get_env();
  problems = [10 40 50 60 70];
  problems = problems(5)
  srcs = eye(3); rcvs = eye(3);
  i_max = 4;
  k = 1;
  for (ip = 1:numel(problems))
    for (is = 1:size(srcs, 1))
      for (ir = 1:size(rcvs, 1))
        o = cts_setopts();
        o.final = 1;
        o.bfn = e.save_dir;
        o.problem = problems(ip);
        o.make_base_fn = @(o) cts_mc_make_base_fn(o, srcs(is,:), rcvs(ir,:));
        ss{k} = cts_mc_run(o, i_max, srcs(is,:), rcvs(ir,:));
        k = k + 1;
      end
    end
  end
  pr('----------------------------\n');
  for (k = 1:numel(ss))
    for (j = 1:numel(ss{k}.cmds)) pr('%s', ss{k}.cmds{j}); end
  end
end

function ss = cts_all_a1_plot (ss)
  for (k = 1:numel(ss))
    if (~isfield(ss{k}.rn{1}, 'tau_interp'))
      ss{k} = cts_a1_init(ss{k});
    end
    clf; cts_a1_p1(ss{k}); drawnow;
    suffix = sprintf('_prob%d_s%dr%d',ss{k}.o.problem,...
                     find(ss{k}.ru{1}.p.src_disl) - 1,...
                     find(ss{k}.ru{1}.p.rcv_traction) - 1);
    pr([ss{k}.o.bfn suffix '.eps\n']);
    print('-depsc', ['~/tmp/cts' suffix '.eps']);
  end
end

function s = cts_run_final (testno)
% Collect together the code for the final tests.
  e = get_env();
  o = cts_setopts();
  o.final = 1;
  o.bfn = e.save_dir;
  i_max = 4;
  switch (testno)
    case 1
      % Periodic (fullspace).
      o.problem = 10;
    case 2
      % W-E periodic, N-S v-BC.
      o.problem = 40;
    case 3
      % W-E periodic, S v-BC, N f.s.
      o.problem = 50;
    case 4
      % f.s. N, v-BC on all other sides.
      o.problem = 60;
    case 5
      % v-BC on all sides.
      o.problem = 70;
  end
  s = cts_mc_run(o, i_max, [1 0 0], [1 0 0]);
  save([o.make_base_fn(o) '_s.mat'], 's');
end

function ttl = cts_get_title (s)
  switch (s.o.problem)
    case 10
      ttl = 'W-E, N-S periodic; fullspace';
    case 40
      ttl = 'W-E periodic; N-S velocity BC; halfspace';
    case {50}
      ttl = 'W-E periodic; S velocity BC; N free surface; halfspace';
    case 60
      ttl = 'E, S, W velocity BC; N free surface; halfspace';
    case 70
      ttl = 'all sides velocity BC; halfspace';
  end
  dirs = {'strike' 'dip' 'normal'};
  src = dirs{find(s.ru{1}.p.src_disl)};
  rcv = dirs{find(s.ru{1}.p.rcv_traction)};
  ttl = sprintf('%s; src: %s; rcv: %s', ttl, src, rcv);
end

function o = cts_setopts ()
% This has the options for the nominal full test.
  e = get_env();
  o.bfn = e.save_dir;
  o.height = 2;
  o.refine = 0;
  o.tol = 1e-6; % was 1e-10;
  o.nthreads = GetNcpu();
  o.neighborhood = 2;
  o.use_bmb_lra = 1;
  o.interp_type = 3;
  o.do_uniform = 0;
  o.problem = 50;
  o.make_base_fn = @cts_make_base_fn;
  o.final = 0;
end

function s = cts_run (o, i_max)
  p = cts_make_problem(o);
  
  base_nbrhd = 1;
  s.o = o; s.rn = {}; s.ru = {}; cmds = {};
  for (i = 0:i_max)
    o.refine = i;
    o.neighborhood = (2^i*(2*base_nbrhd + 1) - 1)/2;

    if (i < i_max)
      o.do_uniform = 0;
      [r c] = cts_write_kvfs(o, p, 1);
      cmds{end+1} = c{3};
      for (j = 1:2) mysystem(c{j}); end
      rects = dc3dm.mRects(r.cb.build_write_filename);
      md = dc3dm.mData(rects);
      r.nelem = size(rects, 2);
      pr('# elem %d\n', r.nelem);
      s.rn{i+1} = r;
    end

    if (i == 0)
      pu = p;
      pu.max_len = diff(p.x([1 end]))/sqrt(r.nelem);
    end
    o.do_uniform = 1;
    [r c] = cts_write_kvfs(o, pu, 1);
    cmds{end+1} = c{3};
    for (j = 1:2) mysystem(c{j}); end
    rects = dc3dm.mRects(r.cb.build_write_filename);
    md = dc3dm.mData(rects);
    r.nelem = size(rects, 2);
    pr('# elem %d\n', r.nelem);
    s.ru{i+1} = r;
  end
  
  for (i = 1:numel(s.rn)) pr('%5d %5d\n', s.rn{i}.nelem, s.ru{i}.nelem); end
  for (i = 1:numel(cmds)) pr(cmds{i}); end
  s.cmds = cmds;
end

function s = cts_mc_run (o, i_max, src, rcv)
% Study multicomponent matrices.
  p = cts_make_problem(o);
  p.src_disl = src;
  p.rcv_traction = rcv;

  base_nbrhd = 1;
  s.o = o; s.rn = {}; s.r0 = {}; s.ru = {}; cmds = {};
  for (i = 0:i_max)
    o.refine = i;
    o.neighborhood = (2^i*(2*base_nbrhd + 1) - 1)/2;
    
    if (i < i_max)
      o.do_uniform = 0;
      [r c] = cts_write_kvfs(o, p, 1);
      cmds{end+1} = c{3};
      for (j = 1:2) mysystem(c{j}); end
      
      cc = r.cc;
      cc.estimate_B_frobenius_only = 1;
      cc.src_disl = [1 0 0];
      cc.rcv_traction = [1 0 0];
      dc3dm.WriteKvf(cc.kvf, cc, true);
      mysystem(c{3});
      com = dc3dm.ReadKvf([cc.hm_write_filename '.compress']);

      r.cc.B_frobenius = com.B_frobenius_estimate;
      dc3dm.WriteKvf(r.cc.kvf, r.cc, true);
            
      rects = dc3dm.mRects(r.cb.build_write_filename);
      md = dc3dm.mData(rects);
      r.nelem = size(rects, 2);
      pr('# elem %d\n', r.nelem);
      s.rn{i+1} = r;
      
      o.neighborhood = 0;
      [r c] = cts_write_kvfs(o, p, 1);
      cmds{end+1} = c{3};
      for (j = 1:2) mysystem(c{j}); end
      r.cc.B_frobenius = com.B_frobenius_estimate;
      dc3dm.WriteKvf(r.cc.kvf, r.cc, true);
      r.nelem = size(rects, 2);
      s.r0{i+1} = r;
    end

    if (i == 0)
      pu = p;
      pu.max_len = diff(p.x([1 end]))/sqrt(r.nelem);
    end
    o.do_uniform = 1;
    o.neighborhood = 0;
    [r c] = cts_write_kvfs(o, pu, 1);
    r.cc.B_frobenius = com.B_frobenius_estimate;
    dc3dm.WriteKvf(r.cc.kvf, r.cc, true);
    cmds{end+1} = c{3};
    for (j = 1:2) mysystem(c{j}); end
    rects = dc3dm.mRects(r.cb.build_write_filename);
    md = dc3dm.mData(rects);
    r.nelem = size(rects, 2);
    pr('# elem %d\n', r.nelem);
    s.ru{i+1} = r;
  end
  
  for (i = 1:numel(s.rn)) pr('%5d %5d\n', s.rn{i}.nelem, s.ru{i}.nelem); end
  for (i = 1:numel(cmds)) pr(cmds{i}); end
  s.cmds = cmds;
end

function s = cts_a1_init (s, varargin)
  o = popt(varargin, {{'wrt' numel(s.ru)}});
  s.ru{o.wrt} = cts_a1_init_r(s, s.ru{o.wrt}, o.wrt, 1);
  for (i = 1:min(numel(s.rn), o.wrt))
    pr('%d ', i);
    s.rn{i} = cts_a1_init_r(s, s.rn{i}, o.wrt);
    s.ru{i} = cts_a1_init_r(s, s.ru{i}, o.wrt);
    if (isfield(s, 'r0'))
      s.r0{i} = cts_a1_init_r(s, s.r0{i}, o.wrt);
    end
  end
  pr('\n');
end

function r = cts_a1_init_r (s, r, wrt, are_same)
  if (nargin < 4) are_same = 0; end
  r = cts_a1_mvp(r, s.ru{wrt}, are_same);
  d = dir([r.cc.hm_write_filename '.hm']);
  r.hm_bytes = d.bytes;
end

function cts_a1_p1 (s)
  function re = relerr1 (a, b)
    n = size(a, 1);
    is = round(n/3):round(2*n/3);
    re = norm(a(is, is) - b(is, is), 'fro') / norm(a(is, is), 'fro');
  end
  
  assert(isfield(s.rn{1}, 'tau_interp'));
  for (i = 1:numel(s.rn))
    ns(1,i) = s.ru{i}.nelem;
    ns(2,i) = s.rn{i}.nelem;
    re(1,i) = relerr(s.ru{end}.tau_interp, s.ru{i}.tau_interp);
    re(2,i) = relerr(s.ru{end}.tau_interp, s.rn{i}.tau_interp);
    if (isfield(s, 'r0'))
      ns(3,i) = s.r0{i}.nelem;
      re(3,i) = relerr(s.ru{end}.tau_interp, s.r0{i}.tau_interp);
    end
  end
  
  ha = 'HorizontalAlignment';
  va = 'VerticalAlignment';
  vec = @(x) x(:);

  function h = contour_wrap (I, N)
    [~, h] = contour(I, N); hold all;
    axis equal;
  end
  function h = contour_signed (I, N)
    hi = imagesc(I); hold all;
    [~, h0] = contour(I >= 0, 1, 'w-');
    [~, hg] = contour(I.*(I > 0), 3, 'k-'); hold off;
    set(hg, 'color', 0.5*[1 1 1]);
    h = [hi(:); h0; hg];
    axis equal;
  end
  imfn = @(I, N) imagesc(I);
  imfn_signed = @contour_signed;
  cmap = gray(1024);
  cmap = cmap(end:-1:1,:);
  colormap(cmap);
  clr = 'kkkk';

  label_letter_wrap = @(varargin) [];
  
  dirs = {'Along-strike shear' 'Along-dip shear' 'Normal'};
  src = dirs{find(s.ru{1}.p.src_disl)};
  rcv = dirs{find(s.ru{1}.p.rcv_traction)};
  
  hs = [];
  pad = 0.1;
  len = (1 - 3*pad)/2;
  fn = cts_a1_make_function(s.ru{end}.p);
  axes('position', [pad, 0.5*(1 + 0*pad), len, len]);
  m = s.ru{end}.m;
  ht = imfn(fn(m.X, m.Y), 11);
  %hs = [hs; ht];
  axis xy; axis equal;
  xlim([0 size(s.ru{end}.tau_interp,2)]);
  ylim([0 size(s.ru{end}.tau_interp,1)]);
  set(gca, 'xtick', [], 'ytick', []);
  %hs = [hs; vec(mytitle(sprintf('Test %s slip', lower(src))))];
  hs = [hs; vec(xlabel(sprintf('Test %s slip', lower(src))))];
  hs = [hs; label_letter_wrap('a')];

  axes('position', [0.5*(1 + pad), 0.5*(1 + 0*pad), len, len]);
  ht = imfn_signed(s.ru{end}.tau_interp, 5);
  %hs = [hs; ht];
  set(gca, 'xtick', [], 'ytick', []);
  axis xy; at; axis equal;
  xlim([0 size(s.ru{end}.tau_interp,2)]);
  ylim([0 size(s.ru{end}.tau_interp,1)]);
  %hs = [hs; vec(mytitle(sprintf('Computed\n%s traction', lower(rcv))))];
  hs = [hs; vec(xlabel(sprintf('Computed %s traction', lower(rcv))))];
  hs = [hs; label_letter_wrap('b')];

  axes('position', [pad, pad, len, len]);
  pat = {'o-' 'v-' 's-'};
  % 2 because we're on a 2D manifold.
  ooa = -2*diff(log(re'))./diff(log(ns'));
  for (i = 1:size(re,1))
    hs = [hs; vec(plot(log10(ns(i,:)'), log10(re(i,:)'), [clr(i) pat{i}]))];
    hold all;
  end
  grid on;
  hs = [hs; vec(xlabel('log_{10} Number of elements'))];
  hs = [hs; vec(mytitle('log_{10} Relative error'))];
  hs = [hs; label_letter_wrap('c')];
  set(gca, 'yaxislocation', 'right');
  xlim([2 5]);

  axes('position', [0.5*(1 + pad), pad, len, len]);
  for (i = 1:size(re,1))
    hs = [hs; vec(plot(log10(ns(i,1:end-1)'), ooa(:,i), [clr(i) pat{i}]))];
    hold all;
  end
  grid on;
  hs = [hs; vec(xlabel('log_{10} Number of elements'))];
  hs = [hs; vec(mytitle('Empirical order of accuracy'))];
  hs = [hs; label_letter_wrap('d')];
  lgd = {'uniform', 'AIGA'};
  if (size(re, 1) == 3) lgd{3} = 'DDMu(n)'; end
  legend(lgd, 'location', 'west');
  ylim([0 2.5]);
  if (0)
    axes('position', [0 0 1 1], 'visible', 'off');
    pr('!!! no title\n');
    ht = text(0.5, 0.5, cts_get_title(s));
    set(ht, 'FontSize', 12, ha, 'center', va, 'middle');
  end
  
  make_pretty(hs);
end

function cts_a1_show_spiky (s)
% Make a figure that shows dramatically (since it is dramatic) the difference
% between DDMu(n) and IGA.
  n = numel(s.r0);
  
  p = get(gcf, 'position');
  p(3:4) = min(p(3:4))*[1 1];
  set(gcf, 'position', p);

  hs = zeros(n, 4);
  bpad = 0.05;
  function spfn (row, col)
    ygap = 0.01;
    dy = (1 - bpad - n*ygap)/n;
    xgap = 0.01;
    dx = (1 - 5*xgap)/4;
    p = [xgap + (col - 1)*(dx + xgap), 1 - row*(dy + ygap), dx, dy];
    hs(row, col) = axes('position', p);
  end
  
  function imfn (img)
    imagesc(s.rn{end}.m.ux, s.rn{end}.m.uy, img);
    axis xy; set(gca, 'xtick', [], 'ytick', []);
    yl = ylim(); xl = xlim();
    if (0)
      % This was for the old problem 60, 70
      %   fn = @(x, y) exp(-20.*(x.^2 + y.^2)).
      ylim(yl(1) + [0.25 0.75]*diff(yl));
      xlim(xl(1) + [0.25 1]*diff(xl));
    else
      % Now I'm using the very tapered crack with a large R.
      ylim(yl(1) + [0 0.75]*diff(yl));
      xlim(xl(1) + [0.25 1]*diff(xl));
    end
    axis equal;
  end
  
  function imlfn (img)
    imfn(log10(abs(img)));
  end

  den = max(abs(s.ru{end}.tau_interp(:)));
  %den = s.ru{end}.tau_interp;
  for (i = 1:n)
    spfn(i, 1); imfn(s.r0{i}.tau_interp);
    spfn(i, 2); imlfn((s.r0{i}.tau_interp - s.ru{end}.tau_interp)./den);
    spfn(i, 3); imfn(s.rn{i}.tau_interp);
    spfn(i, 4); imlfn((s.rn{i}.tau_interp - s.ru{end}.tau_interp)./den);
  end
  
  hss = {hs(:,[1 3]) hs(:,[2 4])};
  for (i = 1:2)
    ca = [inf -inf];
    for (j = 1:numel(hss{i}))
      axes(hss{i}(j));
      ca1 = caxis();
      ca = [min(ca(1), ca1(1)), max(ca(2), ca1(2))];
    end
    if (i == 2)
      ca
      ca = [-5 -1];
    end
    for (j = 1:numel(hss{i}))
      axes(hss{i}(j));
      caxis(ca);
    end
  end
  
  use_gray = 1;
  
  axes('position', [0 0 1 1], 'Visible', 'off');
  top = 0.8*bpad;
  if (use_gray)
    er = sprintf('%d (light) to %d (dark)', ca(1), ca(2));
  else
    er = sprintf('%d to %d', ca(1), ca(2));
  end
  %text(3/8, top, sprintf('log_{10} Error, %s', er));
  ht = [text(1/8, top, 'Traction (DDMu(n))');
        text(3/8, top, 'log_{10} Error');
        text(5/8, top, 'Traction (AIGA)');
        text(7/8, top, 'log_{10} Error')];
  set(ht, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
          'FontSize', 12);
  if (use_gray)
    cmap = gray(1024);
    cmap = cmap(end:-1:1,:);
    colormap(cmap)
  else
    colormap(jet(1024));
  end
  
  if (1)
    drawnow;
    axes(hs(1,1));
    e = explain('make_e', s.rn{1}.cb.build_write_filename);
    ol = explain('ol_make', e);
    h = explain('ol_draw', ol);
    set(h, 'color', 'r');
  end
end

function cts_a1_plot_nvsmem (rs)
  for (i = 1:numel(rs))
    x(i) = rs{i}.nelem;
    y(i) = rs{i}.hm_bytes;
  end
  plot(log10(x), log10(y), '.-',...
       log10(x), log10(x/x(1)*y(1)), '--',...
       log10(x), log10(x.^2/x(1)^2*y(1)), '--');
end

function cts_a1_show_r (r, r_ref)
  if (nargin < 2) r_ref = []; end
  if (isempty(r_ref)) nsp = 2; else nsp = 3; end
  imfn = @(I) imagesc(r.m.ux, r.m.uy, I);
  h(1) = sp(nsp,2,1); imfn(r.tau); cb; ca = caxis(); axis xy;
  h(2) = sp(nsp,2,2); imfn(r.tau_downterp); cb; caxis(ca); axis xy;
  h(3) = sp(nsp,2,3); imfn(log10(abs(r.tau - r.tau_downterp))); cb; axis xy;
  if (~isempty(r_ref))
    imfn = @(I) imagesc(r.m_up.ux, r.m_up.uy, I);
    h(4) = sp(nsp,2,4); imfn(log10(abs(r.tau_interp - r_ref.tau))); cb; axis xy;
    h(5) = sp(nsp,2,5); imfn(r.tau_interp); cb; ca = caxis(); axis xy;
    h(6) = sp(nsp,2,6); imfn(r_ref.tau); cb; caxis(ca); axis xy;
  end
  linkaxes(h);
end

function f = cts_a1_show_function (r)
  x = linspace(-1, 1, 100);
  [X Y] = meshgrid(x, x);
  fn = cts_a1_make_function(r.p);
  f = fn(X,Y);
  sp(211); imagesc(x, x, f); axis xy;
  sp(212); plot(x, [f(:,[1 end]) f([1 end],:)']);
end

function cts_test_fs (r, interp_type)
  rects = dc3dm.mRects(r.cb.build_write_filename);
  [x y] = dc3dm.mCC(rects);
  fn = @(x, y) y; % basic test for exactness
  fn = @(x, y) cos(2*pi*x) + sin(1.5*pi*y);
  f = fn(x, y);
  
  x = linspace(-1, 1, 100);
  y = linspace(-1, 1, 200);
  [X Y] = meshgrid(x, y);
  ft = fn(X, Y);
  
  switch (interp_type)
    case 1
      fi = dc3dm.mLinterpWExtrap(r.cb.build_write_filename, f, X, Y);
    case 3
      fi = dc3dm.mCinterpWExtrap(r.cb.build_write_filename, f, X, Y);
    otherwise
      error('nope');
  end
  h(1) = sp(221); imagesc(x, y, fi);
  h(2) = sp(222); imagesc(x, y, log10(abs((ft - fi)/max(abs(ft(:)))))); cb;
  linkaxes(h);
  sp(212); plot(x, ft([1 end],:), 'b-', x, fi([1 end],:), 'r-');
end

function fn = cts_a1_make_function (p)
  switch (p.problem)
    case 10
      fn = @(x, y) cos(pi*x).*sin(pi*y);
    case {20 30 50 80}
      fn = @(x, y)...
           (cos(1*pi*y/diff(p.y([1 end]))).*...
            exp(-10*y.^2)) +...
           (sin(0*pi + 2*pi*x/diff(p.x([1 end]))).*...
            exp(-10*y.^2)) +...
           sigm.Sigmoid(y, 0, 10, -1, 1, 0, 0.1);
    case 40
      fn = @(x, y)...
           (cos(1*pi*y/diff(p.y([1 end]))).*...
            exp(-10*y.^2)) +...
           (sin(2*pi*x/diff(p.x([1 end]))).*...
            exp(-10*y.^2));
    case {60 70}
      fn = @(x, y) tapered_slip_fn(x, y, 0.8);
  end
end

function slip = tapered_slip_fn (x, y, R)
  r = sqrt(x.^2 + y.^2);
  slip = (1 - (r/R).^2).^(3);
  slip(r >= R) = 0;
end

function r = cts_a1_mvp (r, r_mesh, are_same)
% Compute the MVP for r and interp the result to the mesh given by r_mesh.
  assert(are_same || isfield(r_mesh, 'tau'));
  if (nargin < 3) are_same = 0; end
  
  r.m = make_mesh(r);
  r.m_up = make_mesh(r_mesh);
  
  fn = cts_a1_make_function(r.p);
  if (1)
    rects = dc3dm.mRects(r.cb.build_write_filename);
    [x y] = dc3dm.mCC(rects);
    slip = fn(x(:), y(:));
    md = dc3dm.mData(rects);
    bdy_x = [md.xlim(2) mean(md.xlim) md.xlim(1) mean(md.xlim)];
    bdy_y = [mean(md.ylim) md.ylim(2) mean(md.ylim) md.ylim(1)];
    bdy_vals = fn(bdy_x, bdy_y);
  else
    [slip bdy_vals] = avg_fn_over_elems(fn, r.m, r.m_up);
  end
  bdy_vals(:) = 0;
  
  id = hmmvp('init', [r.cc.hm_write_filename '.hm']);
  r.tau = hmmvp('mvp', id, slip);
  hmmvp('cleanup', id);
  
  bc = dc3dm.ReadBoundaryConditions(r.cc.hm_write_filename);
  r.tau = r.tau + bc*bdy_vals(:);
  
  if (~are_same)
    r.tau_interp = reshape(dc3dm.mCinterpWExtrap(...
      r.cb.build_write_filename, r.tau, r.m_up.X, r.m_up.Y), size(r.m_up.X));
    r.tau_downterp = dc3dm.mCinterpWExtrap(...
      r_mesh.cb.build_write_filename, r_mesh.tau_nop, r.m.x, r.m.y);
    r.tau_downterp = r.tau_downterp(r.m.p);
  else
    r.tau_interp = reshape(r.tau(r.m_up.p), size(r.m_up.X));
    r.tau_downterp = reshape(r.tau(r.m_up.p), size(r.m_up.X));
  end
  
  r.tau_nop = r.tau;
  r.tau = r.tau(r.m.p);
end

function r = test_cts_a1_mvp (r, r_mesh, are_same)
  r = cts_a1_mvp(r, r_mesh, are_same);
  cts_a1_show_function(r); sp(211); cb; ylabel('slip');
  sp(212); imagesc(r.m.ux, r.m.uy, r.tau_interp); ylabel('tau'); cb;
end

function [v bdy_v] = avg_fn_over_elems (fn, m, m_up)
  ids = dc3dm.mIds(m.rm_fn, m_up.X(:), m_up.Y(:));  
  f = fn(m_up.X(:), m_up.Y(:));
  v = zeros(numel(m.x), 1);
  ns = v;
  for (i = 1:numel(f))
    ns(ids(i)) = ns(ids(i)) + 1;
    v(ids(i)) = v(ids(i)) + f(i);
  end
  v = v./ns;
  
  bdy_x = [m.md.xlim(2) mean(m.md.xlim) m.md.xlim(1) mean(m.md.xlim)];
  bdy_y = [mean(m.md.ylim) m.md.ylim(2) mean(m.md.ylim) m.md.ylim(1)];
  bdy_v = fn(bdy_x, bdy_y);
end

function bfn = cts_make_base_fn (o)
  bfn = sprintf(...
    '%s_hgt%duni%drfn%dtol%dnbr%1.2fin%dfin%dprob%d',...
    o.bfn, o.height, o.do_uniform, o.refine, -round(log10(o.tol)),...
    o.neighborhood, o.interp_type, o.final, o.problem);
end

function bfn = cts_mc_make_base_fn (o, src_disl, rcv_traction)
  bfn = cts_make_base_fn(o);
  bfn = sprintf('%s_s%dr%d', bfn, find(src_disl) - 1, find(rcv_traction) - 1);
end

function p = cts_make_problem (o)
  n = 101;
  p.x = linspace(-1, 1, n);
  p.y = linspace(-1, 1, n);
  [X Y] = meshgrid(p.x, p.y);
  base_len = 0.25;
  p.f = base_len*ones(n);

  for (i = 1:o.height)
    p.f(abs(X + 0.9*Y) <= 0.05*base_len*(1/3)^(i/2)) = base_len/2^i;
    p.f(Y > 0 & abs(Y) <= 0.5*base_len*(1/3)^(i/2)) = base_len/2^i;
  end    

  e1 = [1 0 0];
  p.problem = o.problem;
  switch (o.problem)
    case 10
      p.do_fullspace = 1;
      p.bc.ewpbc = 0;
      p.bc.nspbc = 1;
      p.src_disl = e1; p.rcv_traction = e1;
    case {20 30 50}
      p.depth_min = 0;
      p.dipdeg = 12;
      p.bc.ewpbc = 1;
      p.bc.svbc = 2;
      p.src_disl = e1; p.rcv_traction = e1;
    case {40}
      p.depth_min = 1;
      p.dipdeg = 12;
      p.bc.ewpbc = 1;
      p.bc.svbc = 2;
      p.bc.nvbc = 3;
      p.src_disl = e1; p.rcv_traction = e1;
    case {60}
      p.depth_min = 0;
      p.dipdeg = 12;
      p.bc.wvbc = 1;
      p.bc.svbc = 2;
      p.bc.evbc = 3;
      p.src_disl = e1; p.rcv_traction = e1;
    case {70}
      p.depth_min = 1;
      p.dipdeg = 12;
      p.bc.wvbc = 1;
      p.bc.svbc = 2;
      p.bc.evbc = 3;
      p.bc.nvbc = 4;
      p.src_disl = e1; p.rcv_traction = e1;
    otherwise
      error(sprintf('Not a problem: %d.\n', o.problem));
  end
  if (o.final) p.bc.bc_periodic_nlayers = 3;
  else p.bc.bc_periodic_nlayers = 3; end
  
  p.mu = 1;
  p.nu = 0.25;
end

function [r cmds] = cts_write_kvfs (o, p, silent)
  if (nargin < 3) silent = 0; end
  e = get_env();
  cm = cts_write_mesh_kvf(o, p);
  cb = cts_write_build_kvf(o, p);
  cc = cts_write_compress_kvf(o, p);
  cs = {cm cb cc};
  for (i = 1:numel(cs)) kvf_fns{i} = cs{i}.kvf; end
  % Run each of these three commands on the shell command line.
  for (i = 1:3)
    cmds{i} = sprintf('time nice -n 19 %s %s.kvf;\n', e.bin_dc3dm, kvf_fns{i});
    if (~silent) pr(cmds{i}); end
  end
  r.p = p; r.cm = cm; r.cb = cb; r.cc = cc;
end

function c = cts_write_mesh_kvf (o, p)
  c.mesh_write_filename = o.make_base_fn(o);
  if (isfield(p, 'min_len')) c.min_len = p.min_len; else c.min_len = 0; end
  if (isfield(p, 'max_len'))
    c.max_len = p.max_len;
  else
    c.max_len = min(diff(p.x([1 end])), diff(p.y([1 end])))/8;
  end
  c.x = p.x;
  c.y = p.y;
  c.f = p.f;
  if (isfield(o, 'do_uniform')) c.do_uniform = o.do_uniform; end
  c.command = 'mesh';
  c.kvf = [c.mesh_write_filename '_m'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

function c = cts_write_build_kvf (o, p)
  bfn = o.make_base_fn(o);
  c.mesh_read_filename = bfn;
  c.build_write_filename = bfn;
  if (isfield(o, 'refine')) c.refine = o.refine; end
  c = transfer_fields(c, p, {'dipdeg' 'depth_min' 'do_fullspace'});
  c.neighborhood = o.neighborhood;
  c = transfer_fields(c, p.bc, fieldnames(p.bc));
  c.command = 'build';
  c.kvf = [bfn '_b'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

function c = cts_write_compress_kvf (o, p, cb)
  bfn = o.make_base_fn(o);
  c = transfer_fields([], p, {'mu' 'nu' 'src_disl' 'rcv_traction'});
  c.build_read_filename = bfn;
  c.tol = o.tol;
  c.hm_write_filename = bfn;
  c.allow_overwrite = 1;
  c.nthreads = o.nthreads;
  c.command = 'compress';
  c.kvf = [bfn '_c'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

% ------------------------------------------------------------------------------
% Utils.

function mysystem (varargin)
  cmd = sprintf(varargin{:});
  %pr([cmd '\n']);
  %[s r] =...
      system(cmd);
end

function m = make_mesh (r)
% mesh = mesh_init (r)
%   Establish the element order and mesh data.
%   mesh.x, mesh.y are element centers in 'H-matrix order'.
%   mesh.ux, mesh.uy are the unique x and y coordinates.
%   mesh.X, mesh.Y are in mesgrid ordering, which is convenient for making
% images. Call this 'user order'.
%   mesh.p is a permutation that maps H-matrix order to user order. Use it to
% go from matrix-vector product quantities to the order suitable for an image.
  m.rm_fn = r.cm.mesh_write_filename;
  rs = dc3dm.mRects(m.rm_fn);
  [m.x m.y] = dc3dm.mCC(rs);
  m.md = dc3dm.mData(rs);
  m.xlim = m.md.xlim;
  m.ylim = m.md.ylim;

  m.ux = CC(m.md.xlim(1) : m.md.dx : m.md.xlim(2));
  m.uy = CC(m.md.ylim(1) : m.md.dy : m.md.ylim(2));
  [m.X m.Y] = meshgrid(m.ux, m.uy);
  m.p = dc3dm.mIds(r.cm.mesh_write_filename, m.X, m.Y);
end

function make_pretty (h)
  h(end+1) = gca;
  h(end+1) = gcf;
  for (i = 1:length(h))
    try
      if (isfield(get(h(i)), 'LineWidth'))
        set(h(i),'LineWidth', 2, 'MarkerSize', 8);
      end
      if (isfield(get(h(i)), 'FontWeight'))
        set(h(i), 'FontSize', 12, 'fontname', 'helvetica');
      end
    catch, end
  end
end

function ht = mytitle (ttl)
  yl = ylim();
  if (0)
    ht = text(mean(xlim()), yl(2) - 0.03*diff(yl), ttl);
    set(ht, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
  else
    ht = text(mean(xlim()), yl(1) + 0.0*diff(yl), ttl);
    set(ht, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
  end
  set(gca, 'ylimmode', 'manual');
end

function h = label_letter(lt)
  yl = ylim();
  xl = xlim();
  h = text(xl(1), yl(2), ['(' lt ')']);
  set(h, 'horizontalalignment', 'left', 'verticalalignment', 'top');
end

function sd = transfer_fields (sd, ss, flds)
  if (~iscell(flds)) flds = {flds}; end
  for (i = 1:numel(flds))
    if (isfield(ss, flds{i})) sd.(flds{i}) = ss.(flds{i}); end
  end
end
function c = CC (v)
% Cell-centered from node-centered.
  c = 0.5*(v(1:end-1) + v(2:end));
end
function varargout = sp (varargin)
  [varargout{1:nargout}] = subplot(varargin{:});
end
function varargout = pr (varargin)
  fprintf(1, varargin{:});
end
function s = rmsp (s)
  s(s == ' ') = [];
end
function b = allsame (x)
  b = all(x(:) == x(1));
end
function axt (a, brdr)
  if (nargin < 2) brdr = .05; end
  if (nargin < 1) a = gca; end
  axis tight;
  xlim = get(a, 'xlim');
  d = brdr*diff(xlim);
  xlim = [xlim(1)-d, xlim(2)+d];
  set(a, 'xlim', xlim);
  ylim = get(a, 'ylim');
  d = brdr*diff(ylim);
  ylim = [ylim(1)-d, ylim(2)+d];
  set(a, 'ylim', ylim);
end
function varargout = len (varargin)
  [varargout{1:nargout}] = numel(varargin{:});
end
