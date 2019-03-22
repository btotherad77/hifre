function [ header, sequences ] = a_readAll( input_folder )
%args = 'input_folder'  output = [ header, output ]

sequences=' ';
x=dir(input_folder);
[a,b]=size(x);
[h,s,q]=fastqread(strcat(input_folder,'/',x(3,1).name));
for i = 4:a
    x(i,1).name;
    [h0,s0,q0]=fastqread(strcat(input_folder,'/',x(i,1).name));
    h=[h,h0];
    s=[s,s0];
    q=[q,q0];
end
header=h;
sequences=s;

end

