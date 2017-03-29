classdef dc3dm, methods (Static)
% Matlab interface to dc3dm.
%
% dc3dm.WriteKvf(filename, c, allow_overwrite)
%   c is a key-value struct containing string and numeric fields. Write it to
% the file filename.kvf. The file is not necessarily portable to other
% machines. If it exists already, it won't be overwritten unless allow_overwrite
% is true.
%
% c = dc3dm.ReadKvf('Read', filename)
%   Load a key-value struct from the file filename.
%
% rid = mRead(rm_fn)
%   Read the mesh resulting from 'mesh' or 'build' operations. rid is a
% pointer to the assoicated data structure.
%
% dc3dm.mClear(rid)
%   Clear the memory associated with the pointer rid.
%
% dc3dm.mClearAll()
%   Clear the memory associated with all pointers.
%
% rs = dc3dm.mRects(rid)
%   Get a 4 x Nr array of rectangle elements in H-matrix element
% order. rs(1,:) is the SW corner's x coordinate, rs(2,:) is the SW
% corner's y coordinate, rs(3,:) is the x-direction size, and rs(4,:) is
% the y-direction size.
%
% [x y] = dc3dm.mCC(rs)
%   For rs = mRects(rid), get the element centers (x, y) in H-matrix
% element order.
%
% c = dc3dm.mData(rs)
%   For rs = mRects(rid), get the minimum (c.dx, c.dy) and maximum (c.Dx,
% c.Dy) element sizes and the domain limits (c.xlim, c.ylim).
%
% id = dc3dm.mIds(rid, X, Y)
%   For arrays of coordinates (X, Y), return an array of the same size
% containing base-1 indices into the elements.
%
% A = dc3dm.mConstInterp(rid, a, X, Y)
%   'a' is an array of values, one for each element (hence numel(a) ==
% size(rs, 2)). Interpolate by nearest neighbors to the array of coordinates
% (X, Y).
%
% A = dc3dm.mLinterp(rid, a, bdy_vals, X, Y)
%   Interpolate using linear interpolation. bdy_vals has 4 elements in the
% order (east, north, west, south). It sets the values on the boundaries of
% the rectangular fault. These are used only on edges having the velocity
% boundary condition.
%
% A = dc3dm.mCinterp(rid, a, bdy_vals, X, Y)
%   Interpolate using cubic interpolation. This is much slower than mLinterp.
%
% A = dc3dm.mLinterpWExtrap(rid, a, X, Y)
%   Extrapolate at the edges rather than using bdy_vals.
%
% A = dc3dm.mCinterpWExtrap(rid, a, X, Y)
%   Extrapolate at the edges rather than using bdy_vals.
%
% dc3dm.mViewBuild(cb)
%   Show the element sizes as a function of fault position. cb is the
% key-value file Matlab struct for the 'build' operation.
%
% bc = dc3dm.ReadBoundaryConditions(hm_fn)
%   Read the .bc file from the 'compress' operation.

% ------------------------------------------------------------------------------
% Read and write key-value files.
  
  function WriteKvf (fn, c, allow_overwrite)
    if (nargin < 3) allow_overwrite = false; end
    fn = AppendSuffixIfNot(fn, 'kvf');
    kvf('Write', fn, c, allow_overwrite);
  end

  function c = ReadKvf (fn)
  % c = ReadKvf(fn)
  %   Load a struct from the file fn.kvf.
    fptr = fopen(fn, 'rb');
    if (fptr == -1)
      fn = AppendSuffixIfNot(fn, 'kvf');
      fptr = fopen(fn, 'rb');
    end
    if (fptr == -1) error(sprintf('Could not open %s for reading.', fn)); end
    c = kvf('Read', fn);
  end
    
% ------------------------------------------------------------------------------
% Read and analyze mesh, build, and compress data.

  function rid = mRead (rm_fn)
    rid = dc3dm_mex('read', rm_fn);
  end
  
  function mClear (rid)
    c = rid_Init(rid);
    dc3dm_mex('free', c.rid);
    rid_Clear(c);
  end
  
  function mClearAll ()
    dc3dm_mex('freeall');
  end
  
  function rs = mRects (rid)
    c = rid_Init(rid);
    rs = dc3dm_mex('getrects', c.rid);
    rid_Clear(c);
  end
  
  function [x y] = mCC (rs)
    x = rs(1,:) + 0.5*rs(3,:);
    y = rs(2,:) + 0.5*rs(4,:);
    x = x(:);
    y = y(:);
    if (nargout == 1) x = [x y]; end
  end
  
  function c = mData (rs)
    c.dx = min(rs(3,:));
    c.dy = min(rs(4,:));
    c.Dx = max(rs(3,:));
    c.Dy = max(rs(4,:));
    c.xlim = [min(rs(1,:)) max(rs(1,:) + rs(3,:))];
    c.ylim = [min(rs(2,:)) max(rs(2,:) + rs(4,:))];
  end
  
  function id = mIds (rid, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    id = dc3dm_mex('getids', c.rid, X, Y);
    rid_Clear(c);
  end
  
  function A = mConstInterp (rid, a, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    id = dc3dm.mIds(rid, X, Y);
    A = reshape(a(id), size(X));
    rid_Clear(c);
  end
  
  function A = mLinterp (rid, a, bdy_vals, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    A = dc3dm_mex('linterp', c.rid, double(a), double(bdy_vals), X, Y);
    A = reshape(A, size(X));
    rid_Clear(c);
  end
  
  function A = mCinterp (rid, a, bdy_vals, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    A = dc3dm_mex('cinterp', c.rid, double(a), double(bdy_vals), X, Y);
    A = reshape(A, size(X));
    rid_Clear(c);
  end
  
  function A = mLinterpWExtrap (rid, a, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    A = dc3dm_mex('linterpe', c.rid, double(a), X, Y);
    A = reshape(A, size(X));
    rid_Clear(c);
  end
  
  function A = mCinterpWExtrap (rid, a, X, Y)
    assert(all(size(X) == size(Y)));
    c = rid_Init(rid);
    A = dc3dm_mex('cinterpe', c.rid, double(a), X, Y);
    A = reshape(A, size(X));
    rid_Clear(c);
  end
  
  function mViewMesh (c)
  % c is the key-value-file struct for dc3dmMesh.
    rid = dc3dm.mRead(c.mesh_write_filename);
    rs = dc3dm.mRects(rid);
    [X Y] = meshgrid(c.x, c.y);
    ids = dc3dm.mIds(rid, X, Y);
    dc3dm.mClear(rid);
    subplot(221);
    if (isfield(c, 'f'))
      f = c.f;
      f(f < c.min_len) = c.min_len;
      f(f > c.max_len) = c.max_len;
      imagesc(f);
      colorbar; title('Resolution function'); axis image; ca1 = caxis();
    else
      ca1 = [inf -inf];
    end
    subplot(222);
    imagesc(reshape(max(rs(3:4,ids(:))), numel(c.y), numel(c.x)));
    colorbar; title('Element size'); axis image; ca2 = caxis();
    for (i = 1:2)
      subplot(2,2,i); caxis([min(ca1(1), ca2(1)) max(ca1(2), ca2(2))]);
    end
    subplot(223);
    imagesc(f - reshape(max(rs(3:4,ids(:))), numel(c.y), numel(c.x)));
    colorbar; title('Resolution function - Element size');
    for (i = 1:3) h(i) = subplot(2,2,i); end
    linkaxes(h);
  end
  
  function es = mViewBuild (c, n)
  % c is the key-value-file struct for dc3dmBuild.
    if (nargin < 2) n = 1000; end
    rid = dc3dm.mRead(c.build_write_filename);
    rs = dc3dm.mRects(rid);
    md = dc3dm.mData(rs);
    x = CC(linspace(md.xlim(1), md.xlim(2), n+1));
    y = CC(linspace(md.ylim(1), md.ylim(2), n+1));
    [X Y] = meshgrid(x, y);
    id = dc3dm.mIds(rid, X, Y);
    dc3dm.mClear(rid);
    es = reshape(rs(3, id), n, n);
    imagesc(x, y, es);
    colorbar; title('Element size'); axis image;
  end
  
  function bc = ReadBoundaryConditions (hm_fn)
  % bc = ReadBoundaryConditions(hm_fn)
  %   Read the boundary conditions produced by dc3dmCompress. hm_fn is
  % hm_write_filename in the dc3dmCompress key-value file. bc has four
  % columns corresponding to the E, N, W, S boundaries, in that order.
    hm_fn = AppendSuffixIfNot(hm_fn, 'bc');
    fid = fopen(hm_fn, 'rb');
    if (fid == -1) error(sprintf('Could not read %s', hm_fn)); end
    bc = fread(fid, inf, 'double');
    fclose(fid);
    if (mod(numel(bc), 4) ~= 0)
      error('bc file does not have four columns.');
    end
    bc = reshape(bc, numel(bc)/4, 4);
  end
  
  function A = mImage (rid, A, x, y)
    c = rid_Init(rid);
    [X Y] = meshgrid(x, y);
    id = dc3dm_mex('getids', c.rid, X(:), Y(:));
    id = reshape(id, size(X));
    A = A(id);
    rid_Clear(c);
  end
  
end, end

% ------------------------------------------------------------------------------
% Private.

function c = CC (v)
% Cell-centered from node-centered.
  c = (v(1:end-1) + v(2:end))/2;
end

function fn = AppendSuffixIfNot (fn, suffix)
  if (numel(fn) < numel(suffix) + 1) return; end
  if (strcmp(fn(end-numel(suffix):end), ['.' suffix])) return; end
  fn = [fn '.' suffix];
end

function c = rid_Init (rid)
  if (isnumeric(rid) && numel(rid) == 1)
    c.rid = rid;
    c.free = 0;
  elseif (ischar(rid))
    c.rid = dc3dm.mRead(rid);
    c.free = 1;
  else
    error('rid is not in a suitable format.')
  end
end

function rid_Clear (c)
  if (c.free) dc3dm.mClear(c.rid); end
end
