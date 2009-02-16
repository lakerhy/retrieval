function[xout] = testaverage (fname)
  x = load(fname);
  x1 = x(:,6:15);
  xout = mean(x1)
