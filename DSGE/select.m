function [f,ix] = select(xx,endo_names,steady_state)

ix=0;
for ii = 1:length(endo_names)
  name = deblank(endo_names(ii,:));
  TF = strcmp(name, xx);
  if TF == 1
   ix=ii;
   break
  end
end
if ix == 0
    error('fatal (s_state) did not find indicated variable')
end
if exist('steady_state','var')
    f=steady_state(ix);
else
    f=[];
end
