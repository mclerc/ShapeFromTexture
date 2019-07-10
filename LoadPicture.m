function b = LoadPicture(Name);
%N = 512;
if strcmp(Name, 'water') | strcmp(Name, 'sdbp1')| strcmp(Name, 'sdbp2')...
      | strcmp(Name, 'golf')| strcmp(Name, 'sunflowers')...
      | strcmp(Name, 'mud'),
N = 512;
elseif strcmp(Name, 'twocos') | strcmp(Name,'concrete_cylinder') ,
  N = 256;
elseif strcmp(Name, 'sdbp1_128') | strcmp(Name,'twocos_128') ... 
      | strcmp(Name,'golf_128')|strcmp(Name,'golf-retouche') ...
      |strcmp(Name,'pcyl7_128')|strcmp(Name,'gril3I50_26'),
  N = 128;
elseif strcmp(Name, 'Todd1') |strcmp(Name,'Todd2')
  N = 1024;
file = [Name '.raw'];
fid = fopen(file,'r'); 
b = fread(fid,[N N]);
b = b';
elseif strcmp(Name, 'pull_outside_3'),
N = 458;M = 397;
file = [Name '.raw'];
fid = fopen(file,'r'); 
b = fread(fid,[M N]);
b = b';
end
