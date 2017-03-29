function varargout = exmb (varargin)
%DIR Example usage of dc3dm and hmmvp, driven by Matlab.
%   View this file and follow directions prefixed with %DIR in the order they
% are written.
%   Additionally, type 'help hmmvp' and 'help dc3dm' to get detailed
% information on all functionality available in the Matlab interface.
%
%DIR In Matlab,
% >> addpath examples; addpath matlab;
  [varargout{1:nargout}] = feval(varargin{:});
end

%DIR To proceed, you need to build bin/dc3dm.

function r = write_kvfs ()
%DIR Set up a problem and write the dc3dm commands to solve it. On the Matlab
% command line, run
% >> r = ex('write_kvfs');
  o = setopts();
  p = make_props(o);
  
  clf;
  subplot(221); imagesc(p.x, p.y, p.a - p.b); title('a - b'); colorbar;
  subplot(222); imagesc(p.x, p.y, p.a./p.b); title('a/b'); colorbar;
  subplot(223); imagesc(p.x, p.y, p.sigma); title('\sigma'); colorbar;
  subplot(224); imagesc(p.x, p.y, p.h_star); title('h^*_b'); colorbar;
  
  r.cm = write_mesh_kvf(o, p);
  r.cb = write_build_kvf(o);
  r.cc = write_compress_kvf(o, r.cb);
  
  %DIR Run each of these three commands, in order, on the shell command line.
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    fprintf('./bin/dc3dm %s.kvf\n', r.(['c' cmds(i)]).kvf);
  end
end

%DIR To proceed, you need to build matlab/hmmvp and matlab/dc3dm_mex. See
% exmvp_omp.cpp and exmvp_mpi.cpp to see how to compute MVP in C++ and
% exmvp_omp.c for C. The .elem file from the 'dc3dm build' operation provides
% the mesh data and element ordering; the Matlab interface dc3dm is just for
% convenience.

function analyze (r)
%DIR Run
% >> ex('analyze', r);
  
  % Convenience routine to see the element sizes:
  clf;
  dc3dm.mViewBuild(r.cb); axis xy;
  
  % Read the mesh produced by 'dc3dm build'.
  rid = dc3dm.mRead(r.cb.build_write_filename);
  % Get the elements.
  rs = dc3dm.mRects(rid);
  % Element centers:
  [cx cy] = dc3dm.mCC(rs);
  dc3dm.mClear(rid);
  
  % Compare these data with those in the .elem file, which is a plain text
  % file formatted as comma separated values. This format is easy to parse,
  % though Matlab happens to have a built-in reader for it, csvread.
  ers = parse_elem_file(r.cb.build_write_filename);
  % Confirm ers contains equivalent information to rs.
  fprintf('Should be <= ~%1.1e: %1.1e %1.1e %1.1e %1.1e\n', ...
          eps(10), relerr(cx(:), ers(1,:)'), relerr(cy(:), ers(2,:)'), ...
          relerr(rs(3,:), ers(3,:)), relerr(rs(4,:), ers(4,:)));
end

%DIR The next two routines demonstrate MVP functionality without focusing on
% a specific calculation.

function demo_mvp (r)
%DIR To see how to compute a matrix-vector product (MVP) in Matlab, read the
% following code and run (with b from 'build')
% >> ex('demo_mvp', r);
  addpaths();
  hm_filename = r.cc.hm_write_filename;
  % Load the H-matrix A into memory. id is a pointer to this H-matrix. Use 4
  % threads if they are available.
  id = hmmvp('init', hm_filename, 4);
  % Get the matrix's number of rows.
  m = hmmvp('getm', id);
  % Get the number of columns.
  n = hmmvp('getn', id);
  % Make a test vector.
  x = randn(n, 1);
  % Compute y = A*x.
  y = hmmvp('mvp', id, x);
  % y should be equal to size(A, 1).
  assert(numel(y) == m);
  % Clean up the memory for A. In practice, don't clean up until you're
  % completely done with this H-matrix.
  hmmvp('cleanup', id);
end

function demo_mvp_advanced (r)
%DIR Products of the form y(rs,:) = A(rs,cs)*x(cs,:) can be computed. Type
% >> ex('demo_mvp_advanced', r);
% to run this demo.
  addpaths();
  hm_filename = r.cc.hm_write_filename;
  % Load the H-matrix A into memory. id is a pointer to this H-matrix. Use 4
  % threads if they are available. Additionally, allow up to 5 MVP in one
  % call.
  id = hmmvp('init', hm_filename, 4, 5);
  % Get the matrix's number of rows.
  m = hmmvp('getm', id);
  % Get the number of columns.
  n = hmmvp('getn', id);
 
  % Make a random set of vectors.
  X = randn(n, 5);
  % Do the basic full MVP:
  Y = hmmvp('mvp', id, X);

  % Now use just a subset of the rows:
  rs = 1:2:m;
  Y1 = hmmvp('mvp', id, X, rs);
  fprintf(['Should be the same to ~sqrt(eps) or eps, ', ...
           'depending on tolerance: %1.1e\n'],...
	  relerr(Y(rs,:), Y1(rs,:)));
  
  % Subset of cols:
  cs = 1:3:n;
  ncs = setdiff(1:n, cs);
  Y2a = hmmvp('mvp', id, X, [], cs);
  Y2b = hmmvp('mvp', id, X, [], ncs);
  fprintf('relerr: %1.1e\n', relerr(Y, Y2a + Y2b));
  
  % Subset of rows and cols.
  Y2a = hmmvp('mvp', id, X, rs, cs);
  Y2b = hmmvp('mvp', id, X, rs, ncs);
  fprintf('relerr: %1.1e\n', relerr(Y(rs,:), Y2a(rs,:) + Y2b(rs,:)));
  
  % For repeated calls to 'mvp' having the same row and column subsets, call
  %     hmmvp('ststate',id);
  % to save the internal state associated with the index sets. For small
  % index sets, this can offer quite a speedup. Call
  %     hmmvp('strelease',id);
  % when finished. NB: There is no error checking, so failing to release the
  % internal state will be a silent bug.
  hmmvp('ststate',id);
  for (i = 1:10)
    Y2a = hmmvp('mvp', id, X, rs, cs);
  end
  hmmvp('strelease',id);
  
  hmmvp('cleanup', id);
end

function r = recompress (r)
%DIR Recompress an H-matrix to increase its compression efficiency. Run
% >> ex('recompress', r);
  addpaths();
  
  % Copy the old key-value struct to a new one, add a use_hmat_filename
  % field, and modify the output filenames.
  r.ccr = r.cc;
  r.ccr.hm_use_filename = r.ccr.hm_write_filename;
  r.ccr.hm_write_filename = [r.ccr.hm_write_filename '_rc.hm'];
  r.ccr.kvf = [r.ccr.kvf '_rc.kvf'];
  % Now make the tolerance looser. It can also be left the same. On larger
  % matrices, recompressing at the same tolerance can slightly improve
  % compression. For any tol >= the old one, recompression should be a lot
  % faster than the original compression.
  r.ccr.tol = 1e2*r.ccr.tol;
  % Write the new kvf-file.
  dc3dm.WriteKvf(r.ccr.kvf, r.ccr, 1);
  
  %DIR In a shell, run this:
  fprintf('./bin/dc3dm %s\n', r.ccr.kvf);
end

function demo_mvp_slip (r)
%DIR Carry out typical operations for a real problem. Run
% >> ex('demo_mvp_slip', r);
  
  % Get element centers. (Could use the .elem file.)
  rid = dc3dm.mRead(r.cb.build_write_filename);
  rs = dc3dm.mRects(rid);
  [cx cy] = dc3dm.mCC(rs);
  md = dc3dm.mData(rs);
  
  % Make a slip distribution. Respect the BCs.
  slip_fn = @(x, y) cos(2*pi*(x + 0.39*diff(md.xlim))/diff(md.xlim)).* ...
            sin(2*pi*y/diff(md.ylim)) + ...
            0.*y/diff(md.ylim);
  slip = slip_fn(cx, cy);
  
  % Get boundary condition data.
  bc = dc3dm.ReadBoundaryConditions(r.cc.hm_write_filename);
  % Boundary values are ordered (east, north, west, south); only those for
  % non-0 velocity-BC matter. Here I fill in all values even though only the
  % S one is necessary.
  bdy_vals = [slip_fn(md.xlim(2), 0), slip_fn(0, md.ylim(2)), ...
              slip_fn(md.xlim(1), 0), slip_fn(0, md.ylim(1))];
  
  % Compute traction.
  id = hmmvp('init', r.cc.hm_write_filename, 4, 1);
  traction = hmmvp('mvp', id, slip) + bc*bdy_vals(:);
  hmmvp('cleanup', id);
  
  % Plot. Show the raw (const interp) mesh values and two different
  % interpolations.
  x = CC(linspace(md.xlim(1), md.xlim(2), round(diff(md.xlim)/md.dx) + 1));
  y = CC(linspace(md.ylim(1), md.ylim(2), round(diff(md.ylim)/md.dy) + 1));
  [X Y] = meshgrid(x, y);
  slip_c = dc3dm.mConstInterp(rid, slip, X, Y);
  traction_c = 1e-6*dc3dm.mConstInterp(rid, traction, X, Y);
  slip_1 = dc3dm.mLinterp(rid, slip, bdy_vals, X, Y);
  traction_1 = 1e-6*dc3dm.mLinterpWExtrap(rid, traction, X, Y);
  slip_3 = dc3dm.mCinterp(rid, slip, bdy_vals, X, Y);
  traction_3 = 1e-6*dc3dm.mCinterpWExtrap(rid, traction, X, Y);
  
  dc3dm.mClear(rid);
  
  clf;
  img = @(im) imagesc(x, y, im);
  h(1) = subplot(321); img(slip_c); title('slip, mesh resolution');
  h(2) = subplot(322); img(traction_c); title('traction, mesh resolution');
  h(3) = subplot(323); img(slip_1); title('slip, linear interp');
  h(4) = subplot(324); img(traction_1); title('traction, linear interp');
  h(5) = subplot(325); img(slip_3); title('slip, cubic interp');
  h(6) = subplot(326); img(traction_3); title('traction, cubic interp');
  linkaxes(h); zoom on;
  for (i = 1:numel(h))
    subplot(h(i));
    if (mod(i, 2) == 0) caxis(135*[-1 1]); else caxis([-1 1]); end
    colorbar;
    axis equal; axis xy; axis tight;
    % Draw the mesh.
    draw_rects_r(rs, 0, 'k');
  end
end

function r = run_test ()
%DIR Run
% >> r = ex('run_test');
% This test runs dc3dm on 3 problems to illustrate its purpose.
  o = setopts();
  o.want_free_surface = 0;
  r.p = make_props(o);

  % AIGA-8: approximate IGA with neighborhood = 8.
  o.neighborhood = 8;
  r.n8.o = o;
  r.n8.cm = write_mesh_kvf(o, r.p);
  r.n8.cb = write_build_kvf(o);
  r.n8.cc = write_compress_kvf(o, r.n8.cb);
  % Run dc3dm.
  cmds = 'mbc';
  for (i = 1:numel(cmds))
    system(sprintf('./bin/dc3dm %s.kvf\n', r.n8.(['c' cmds(i)]).kvf));
  end
  
  % AIGA-0: approximate IGA with neighborhood = 0. This is equivalent to
  % applying the DDMu method on a nonuniform mesh.
  o.neighborhood = 0;
  r.n0.o = o;
  r.n0.cm = write_mesh_kvf(o, r.p);
  r.n0.cb = write_build_kvf(o);
  r.n0.cc = write_compress_kvf(o, r.n0.cb);
  for (i = 1:numel(cmds))
    system(sprintf('./bin/dc3dm %s.kvf\n', r.n0.(['c' cmds(i)]).kvf));
  end
  
  % DDMu on the fine-resolution mesh.
  rid = dc3dm.mRead(r.n8.cb.build_write_filename);
  rs = dc3dm.mRects(rid);
  md = dc3dm.mData(rs);
  dc3dm.mClear(rid);
  o.do_uniform = 1;
  o.neighborhood = 0;
  % Make the element size the smallest in the nonuniform mesh.
  o.max_len = min(diff(md.xlim)/md.dx, diff(md.ylim)/md.dy);
  r.u.o = o;
  r.u.cm = write_mesh_kvf(o, r.p);
  r.u.cb = write_build_kvf(o);
  r.u.cc = write_compress_kvf(o, r.u.cb);
  for (i = 1:numel(cmds))
    system(sprintf('./bin/dc3dm %s.kvf\n', r.u.(['c' cmds(i)]).kvf));
  end
end

function show_test (r)
%DIR Run
% >> ex('show_test', r);
% where r is from 'run_test'. Shows the behavior of the three methods.
    
  % Show results.
  flds = {'u' 'n0' 'n8'};
  for (i = 1:3)
    cb = r.(flds{i}).cb;
    cc = r.(flds{i}).cc;
  
    bc = dc3dm.ReadBoundaryConditions(cc.hm_write_filename);
    
    % Get element centers.
    rid = dc3dm.mRead(cb.build_write_filename);
    rs = dc3dm.mRects(rid);
    [cx cy] = dc3dm.mCC(rs);
  
    if (i == 1)
      md = dc3dm.mData(rs);

      % One mesh for all three.
      x = CC(linspace(md.xlim(1), md.xlim(2), round(diff(md.xlim)/md.dx) + 1));
      y = CC(linspace(md.ylim(1), md.ylim(2), round(diff(md.ylim)/md.dy) + 1));
      [X Y] = meshgrid(x, y);
      
      img = @(im) imagesc(x, y, im);
      
      % The slip function should go smoothly to 0 at N and S boundaries to
      % match the boundary conditions.
      slip_fn = @(x, y) cos(2*pi*(x + 0.1*diff(md.xlim))/diff(md.xlim)).* ...
                exp(-40*((y + 0.05*diff(md.ylim))/diff(md.ylim)).^2);
      
      bdy_vals = [slip_fn(md.xlim(2), 0), slip_fn(0, md.ylim(2)), ...
                  slip_fn(md.xlim(1), 0), slip_fn(0, md.ylim(1))];
    end
    
    slip = slip_fn(cx, cy);
    
    id = hmmvp('init', cc.hm_write_filename, 4, 1);
    traction = hmmvp('mvp', id, slip) + bc*bdy_vals(:);
    hmmvp('cleanup', id);
    
    % Images.
    h(i, 1) = subplot(3, 4, 4*(i-1) + 1);
    if (i == 1)
      slip = dc3dm.mConstInterp(rid, slip, X, Y);
    else
      slip = dc3dm.mCinterp(rid, slip, bdy_vals, X, Y);
    end
    img(slip);
    caxis([-1 1]);
    if (i > 1) draw_rects_r(rs, 0, 'k'); end
    title(sprintf('slip_{%s}', flds{i}));
    
    h(i, 2) = subplot(3, 4, 4*(i-1) + 2);
    if (i == 1)
      traction = dc3dm.mConstInterp(rid, traction, X, Y);
    else
      traction = dc3dm.mCinterpWExtrap(rid, traction, X, Y);
    end
    traction = 1e-6*traction;
    img(traction);
    caxis(166*[-1 1]);
    title(sprintf('traction_{%s}', flds{i}));
    
    % Record the uniform-mesh results for comparison with the other two.
    if (i == 1)
      slip_u = slip;
      traction_u = traction;
    end
    
    h(i, 3) = subplot(3, 4, 4*(i-1) + 3);
    img(abs(slip - slip_u));
    title(sprintf('slip_{%s} - slip_u', flds{i}));
    caxis(0.1*[0 1]);

    % The key observation here is that where the mesh is coarse, both methods
    % have considerable error due simply to coarsenss and interpolation
    % error. But within the borders of the refined region, the naive method
    % also has considerable error, whereas AIGA-8 does not.
    h(i, 4) = subplot(3, 4, 4*(i-1) + 4);
    img(abs(traction - traction_u));
    title(sprintf('traction_{%s} - traction_u', flds{i}));
    caxis(38*[0 1]);
    
    for (j = 1:4)
      subplot(h(i, end-4+j));
      axis equal; axis xy; axis tight;
      set(gca, 'xtick', [], 'ytick', []);
    end
    
    dc3dm.mClear(rid);
  end
  linkaxes(h);
end

% ------------------------------------------------------------------------------
% Implementation details.

function addpaths ()
  addpath('matlab/');
end

function o = setopts ()
% Lengths are in [m].
  o.rfac = 1;
  o.len_fac = 1;
  o.vary_fac = 2;
  o.want_free_surface = 1;
  o.tol = 1e-5;
  o.problem = 1;
  o.do_uniform = 0;
  o.neighborhood = 8;
  o.max_len = inf;
  o.nthreads = 4;
  o.dir = './tmp/';
end

function bfn = make_base_fn (o)
% Come up with an informative base filename.
  bfn = sprintf('%sdc3t_rf%1.2flf%1.2fvf%1.2fnbr%du%d', o.dir, o.rfac, ...
                o.len_fac, o.vary_fac, o.neighborhood, o.do_uniform);
end

function p = make_props (o)
% Set the frictional and other properties for the domain.
  dip_len = o.len_fac*1e3;
  strike_len = o.len_fac*1e3;
  n = 1001;

  p.x = linspace(-0.5*strike_len, 0.5*strike_len, n);
  p.y = linspace(-0.5*dip_len, 0.5*dip_len, n);
  [X Y] = meshgrid(p.x, p.y);
  
  radius = 0.1*o.len_fac*dip_len; % radius of disk
  t_width = 0.5*radius; % transition width
  
  r = sqrt(X.^2 + Y.^2) / radius;
  alpha = calc_transition_width(t_width/radius, 0.99);
  w_sigma = calc_sigmoid(r, 1, alpha, 0, max(r(:)), 0, 1);
  w_amb = calc_sigmoid(r, 2, alpha, 0, max(r(:)), 0, 1);
  
  one = ones(size(X));
  bl = 0.01; bs = 0.005;
  amb_vw = -0.005; amb_vs = 0.005;
  d_c = 1e-4;
  sigma_s = 1e6; sigma_l = o.vary_fac*sigma_s;

  p.mu = 3e10;
  p.nu = 0.25;
  p.b = bl*(1 - w_amb) + bs*w_amb;
  p.a = p.b + amb_vw*(1 - w_amb) + amb_vs*w_amb;
  p.d_c = d_c*one;
  p.sigma = sigma_l*(1 - w_sigma) + sigma_s*w_sigma;
  p.h_star = 1.377*p.mu/(1 - p.nu)*p.d_c./(p.sigma.*p.b);
end

function c = write_mesh_kvf (o, p)
% Write the key-value file for dc3dmMesh.
  c.mesh_write_filename = make_base_fn(o);
  c.do_uniform = o.do_uniform;
  % Min and max element lengths don't really matter unless they are used to
  % bound the resolution function f. Here f is well behaved so I set the min
  % to 0 and make sure there are at least 8 elements in each direction of the
  % domain. (max_len is all that matters if the do_uniform option is true.)
  c.min_len = 0;
  if (isinf(o.max_len))
    c.max_len = min(diff(p.x([1 end])), diff(p.y([1 end])))/8;
  else
    c.max_len = o.max_len;
  end
  % Create a tensor mesh ...
  c.x = p.x;
  c.y = p.y;
  % ... on which f, the resolution function, is set. f has the same units as
  % x and y. Make sure there are 5 o.rfac elements in each direction per h*.
  c.f = p.h_star/(o.rfac*5);
  c.command = 'mesh';
  c.kvf = [make_base_fn(o) '_m'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

function c = write_build_kvf (o)
% Write the key-value file for dc3dmBuild.
  bfn = make_base_fn(o);
  c.mesh_read_filename = bfn;
  c.build_write_filename = sprintf('%s_p%d', bfn, o.problem);
  
  % Setting depth_min to 0 makes the 0-depth boundary (N or S depending on
  % the sign of dipdeg) have a free surface boundary condition.
  c.depth_min = 100;
  if (o.want_free_surface) o.depth_min = 0; end
  switch (o.problem)
    case 1 % subduction
      % Positive dip makes the north boundary be at the surface and the south
      % boundary to have the velocity BC.
      c.dipdeg = 12;
      % The fault is periodic along-strike.
      c.ewpbc = 0;
      % This is the velocity boundary condition at depth.
      c.svbc = 2;
    otherwise
      error(sprintf('%d is not a valid problem number.', o.problem));
  end
  
  c.neighborhood = o.neighborhood;
  c.bc_periodic_nlayers = 3;
  
  c.command = 'build';
  c.kvf = [bfn '_b'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

function c = write_compress_kvf (o, cb)
% Write the key-value file for dc3dmCompress.
  switch (o.problem)
    case 1
      v = [1 2 0];
      v = v/norm(v);
      c.src_disl = v;
      c.rcv_traction = c.src_disl;
      c.component = 1;
    otherwise
      error(sprintf('%d is not a valid problem number.', o.problem));
  end

  bfn = make_base_fn(o);
  c.mesh_read_filename = bfn;
  c.build_read_filename = cb.build_write_filename;
  c.tol = o.tol;
  c.hm_write_filename = sprintf( ...
    '%s_tol%1.1f', cb.build_write_filename, -log10(c.tol));

  c.allow_overwrite = 1;
  c.mu = 3e10;
  c.nu = 0.25;
  c.nthreads = o.nthreads;
  
  c.command = 'compress';
  c.kvf = [bfn '_c'];
  dc3dm.WriteKvf(c.kvf, c, true);
end

function cf = calc_cf (hmfn)
  [id nnz] = hmmvp('init', hmfn);
  m = hmmvp('getm', id);
  n = hmmvp('getn', id);
  hmmvp('cleanup', id);
  cf = m*n/nnz;
end

function sd = transfer_fields (sd, ss, flds)
  if (~iscell(flds)) flds = {flds}; end
  for (i = 1:numel(flds))
    if (isfield(ss, flds{i})) sd.(flds{i}) = ss.(flds{i}); end
  end
end

function re = relerr (a, b, p)
  if (nargin < 3) p = 'fro'; end
  re = norm(a - b, p)/norm(a, p);
  if (isnan(re))
    if (isempty(a) || all(a(:) == b(:))) re = 0; end
  end
end

function y = calc_sigmoid (x, xt, a, xs, xe, ys, ye)
% y = calc_sigmoid(x, xt, alpha, p, xs, xe, ys, ye)
%   xt is the transition point.
%   ys(xs), ye(xe) are the values at reference points.
%   alpha is the constant specifying the transition width. See
% calc_transition_width for more.
  assert(xe > xs);
  fn = @(x) 1./(1 + exp(-a.*x));
  y0s = fn(xs - xt);
  y0e = fn(xe - xt);
  y = ys + (ye - ys).*(fn(x - xt) - y0s)./(y0e - y0s);
end

function alpha = calc_transition_width (width, at_p)
% alpha = calc_transition_width(width, at_p)
%   The transition has width 'width' in the sense that
%       abs(diff(y(x +/- width/2))) = at_p * abs(y1 - ye).
% So at_p should be something like 0.9.
%   We do this calculation for the exponent p = 1 only.
  assert(at_p > 0 && at_p < 1);
  assert(width > 0);
  alpha = -2/width*log(2/(1 + at_p) - 1);
  assert(alpha > 0);
end

function rs = parse_elem_file (fn)
% Parse the .elem file. rs is 4xNr, where Nr is the number of elements in the
% mesh, ordered the same as the H-matrix. Rows 1 and 2 are element center (x, y)
% in the fault coordinates, and rows 3 and 4 are the along-x and along-y
% lengths.
  fid = fopen([fn '.elem'], 'r');
  if (fid == -1) error(sprintf('Can''t read %s.elem', rs)); end
  while (1)
    ln = fgetl(fid);
    if (isempty(ln) || ln(1) < '0' || ln(1) > '9') break; end
    a = sscanf(ln, '%d,%e,%e,%e,%e');
    rs(:, a(1)) = a(2:5)';
  end
  fclose(fid);
end

function c = CC (v)
% Cell-centered from node-centered
  c = 0.5*(v(1:end-1) + v(2:end));
end

function h = draw_rects_r (r, varargin)
  h = [];
  for (i = 1:size(r, 2))
    h1 = draw_rect(r(1,i), r(2,i), r(1,i) + r(3,i), r(2,i) + r(4,i),...
                   varargin{:});
    h = [h; h1(:)];
  end
end

function h = draw_rect (xlo, ylo, xhi, yhi, filled, clr, varargin)
  if (filled)
    h = fill([xlo xhi xhi xlo], [ylo ylo yhi yhi], clr);
  else
    h = line([xlo xhi; xhi xhi; xhi xlo; xlo xlo]',...
             [ylo ylo; ylo yhi; yhi yhi; yhi ylo]');
    set(h, 'color', clr, varargin{:});
  end
end
