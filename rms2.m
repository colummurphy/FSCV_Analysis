function y = rms2(x, varargin)
if isreal(x)
  y = sqrt(mean(x .* x, varargin{:}));
else
  absx = abs(x);
  y    = sqrt(mean(absx .* absx, varargin{:}));
end
end