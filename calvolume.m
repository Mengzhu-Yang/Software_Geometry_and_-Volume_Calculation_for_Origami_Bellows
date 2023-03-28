function v = calvolume(a,b,c,d)
h = [b-a;c-a;d-a];
v = abs(det(h))/6;
end