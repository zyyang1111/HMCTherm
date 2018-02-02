function [ section section_name ] = read_mcpat_power_area( mcpat_output_file )

    fid=fopen(mcpat_output_file,'r');

    %%%Check if file exists
    if fid==-1
        display(['ERROR: ' mcpat_output_file 'is not a valid mcpat output file'])
        return
    end

    flag=0; %flag indicates the preceeding line was "***" which indicates that the next line will be a section header
    count=1;
    line_number=1;
    while ~feof(fid)
        linebuffer=fgetl(fid);
        if ~isempty(regexp(linebuffer,'.*:.*')) || flag==1 %if this line is a section header
            tok_cnt=1;
            remaining=linebuffer;
            while ~isempty(remaining)
                [str{tok_cnt},remaining] = strtok(remaining);
                tok_cnt=tok_cnt+1;
            end
            name=[]; %construct the name of the section (could be multiple tokens)
            for i=1:tok_cnt-1
                if str{i}(end)==':'
                    name=[name str{i}(1:end-1)]; %remove the colon at the end
                    break
                else
                    name=[name str{i} ' '];
                end
            end
            if flag && str{i}(end)~=':'
                name=name(1:end-1); %if there was no colon at the end, remove the final space
            end
            section{count}.section_name=name;
            section_name{count}=name;
            section{count}.names=[];
            section{count}.values=[];
            count=count+1;                
%             display([linebuffer ''])
        elseif count>1 && ~isempty(linebuffer) %otherwise it is a section member
            tok_cnt=1;
            remaining=linebuffer;
            while ~isempty(remaining)
                [str{tok_cnt},remaining] = strtok(remaining);
                tok_cnt=tok_cnt+1;
            end
            name=[];
            for i=1:tok_cnt-1
                if str{i}=='='
                    name=name(1:end-1);
                    break
                else
                    name=[name str{i} ' '];
                end
            end
                if length(str)>=i+1
                    value=str{i+1}; %value equals first token after "="
                else
                    value=[]; %if no "=" or no token after "=" then value is null
                end
            section{count-1}.names=[section{count-1}.names; {name}];
            section{count-1}.values=[section{count-1}.values; {value}];
        end
        line_number=line_number+1;
        
        if ~isempty(linebuffer)
            if linebuffer(1)=='*'; %if line is "****" set flag, else clear it
                flag=1;
%                 linebuffer
            else
                flag=0;
            end
        else
            flag=0;
        end
        
    end

    fclose(fid);

end

