eendtime = 51100;
runtime = 0;
dtout = 10;
dtin = 1;
disp('Initializing chiller minimally...');
chiller(0);
disp('Chiller initialized.');

outfile = fopen('IOFiles\FF_Output.txt','wt');
bcfile  = fopen('IOFiles\FFBCs.txt','rt');
if(outfile==-1)
    disp('Could not open output file for write access.');
    return;
end
if(bcfile==-1)
    disp('Could not open BC file for read access.');
    return;
end

while(runtime<=51100)
    if(runtime==dtout | runtime==0)
        runtime=0;
        [u,count] = fscanf(bcfile,'%f %f %f %f %f',5);
        disp('Boundary conditions refreshed.');
        if(count~=5)
            disp('Error in BC file.');
            clear chiller;
            fclose(bcfile);
            fclose(outfile);
            return;
        end
    end        
    y = chiller(dtin,u);
    y'
    for i=1:29
        fprintf(outfile,'%f\t',y(i));
    end
    fprintf(outfile,'\n');
    runtime = runtime + 1;
end


fclose(bcfile);
fclose(outfile);