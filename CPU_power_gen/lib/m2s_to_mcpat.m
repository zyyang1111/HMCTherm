% function [ mem_file_section, mem_file_section_name, pipeline_file_section, pipeline_file_section_name, memconfig_module, memconfig_module_name, memconfig_network, memconfig_network_name, memconfig_entry, memconfig_entry_name, memconfig_cachegeometry, memconfig_cachegeometry_name ] = m2s_to_mcpat( mem_file, pipeline_file, memconfig_file, net_file, xml_file, num_MC, clk, router_ports, network_latency, NOC_width )
function [  ] = m2s_to_mcpat( mem_file, pipeline_file, memconfig_file, net_file, xml_file, num_MC, clk, router_ports, network_latency, NOC_width, mem_clk, mem_size, num_banks )

    overwrite=1; %if a mcpat.xml already exists, should it be overwritten?

    %x86 constants
    size_arch_FRF=8;
    size_arch_IRF=32;

    %OPTIONS
    L1_dir=0;
    L2_dir=0;
    num_L3=0;
    homo_core=0;
    homo_L3=1;
    homo_L1_dir=1;
    homo_L2_dir=1;
    tech=45;
%     tech=32;
    temp=320;
    icache_geo='geo-i-l1';
    icache_mod='mod-il1-';
    dcache_geo='geo-d-l1';
    dcache_mod='mod-dl1-';   

    for i=1:num_MC
        mem_mod{i}=['mod-mm' num2str(i-1)];
    end
    L2_mod='mod-l2-';
    

    fid_mem=fopen(mem_file,'r');
    fid_pipeline=fopen(pipeline_file,'r');
    fid_network=fopen(net_file,'r');
    fid_memconfig=fopen(memconfig_file,'r');
    fid_xml=fopen(xml_file,'r');

    if fid_mem==-1
        display(['ERROR: ' mem_file ' is not a valid multi2sim memory report'])
        return
    end

    if fid_pipeline==-1
        display(['ERROR: ' pipeline_file ' is not a valid multi2sim pipeline report'])
        return
    end
    if fid_network==-1
        display(['ERROR: ' net_file ' is not a valid multi2sim network report'])
        return
    end
    if fid_memconfig==-1
        display(['ERROR: ' memconfig_file ' is not a valid multi2sim memory configuration'])
        return
    end    
    
    if ~overwrite
        if fid_xml~=-1 %if xml file already exists, just return because the function has been called in the past
            fclose(fid_xml);
            fclose(fid_mem);
            fclose(fid_pipeline);
            fclose(fid_memconfig);
            fclose(fid_network);
            return
        end
    end

    count=1;
    line_number=1;
    while ~feof(fid_mem)
        linebuffer=fgetl(fid_mem);
        if regexp(linebuffer,'\[ .* \]')
            mem_file_section{count}.section_name=linebuffer(3:(length(linebuffer)-2));
            mem_file_section_name{count}=linebuffer(3:(length(linebuffer)-2));
            mem_file_section{count}.names=[];
            mem_file_section{count}.values=[];
            count=count+1;
        elseif count>1
            if ~isempty(linebuffer)
                if linebuffer(1)~=';' %not a comment
                    remaining=linebuffer;
                    tok_cnt=1;
                    while ~isempty(remaining)
                        [str{tok_cnt},remaining] = strtok(remaining);
                        tok_cnt=tok_cnt+1;
                    end
                    for i=1:tok_cnt-1
                        if strcmp(str{i},'=')
                            ptr_name=i-1;
                            ptr_val=i+1;                    
                            mem_file_section{count-1}.names=[mem_file_section{count-1}.names; {str{ptr_name}}];
                            mem_file_section{count-1}.values=[mem_file_section{count-1}.values; {str{ptr_val}}];
                        end
                    end
                end
            end
        end
        line_number=line_number+1;
    end
    
    count=1;
    line_number=1;
    while ~feof(fid_pipeline)
        linebuffer=fgetl(fid_pipeline);
        if regexp(linebuffer,'\[ .* \]')
            pipeline_file_section{count}.section_name=linebuffer(3:(length(linebuffer)-2));
            pipeline_file_section_name{count}=linebuffer(3:(length(linebuffer)-2));
            pipeline_file_section{count}.names=[];
            pipeline_file_section{count}.values=[];
            count=count+1;
        elseif count>1
            if ~isempty(linebuffer)
                if linebuffer(1)~=';'
                    remaining=linebuffer;
                    tok_cnt=1;
                    while ~isempty(remaining)
                        [str{tok_cnt},remaining] = strtok(remaining);
                        tok_cnt=tok_cnt+1;
                    end
                    for i=1:tok_cnt-1
                        if strcmp(str{i},'=')
                            ptr_name=i-1;
                            ptr_val=i+1;                    
                            pipeline_file_section{count-1}.names=[pipeline_file_section{count-1}.names; {str{ptr_name}}];
                            pipeline_file_section{count-1}.values=[pipeline_file_section{count-1}.values; {str{ptr_val}}];
                        end
                    end
                end
            end
        end
        line_number=line_number+1;
    end
    
    count=1;
    line_number=1;
    while ~feof(fid_network)
        linebuffer=fgetl(fid_network);
        if regexp(linebuffer,'\[ .* \]')
            network_file_section{count}.section_name=linebuffer(3:(length(linebuffer)-2));
            network_file_section_name{count}=linebuffer(3:(length(linebuffer)-2));
            network_file_section{count}.names=[];
            network_file_section{count}.values=[];
            count=count+1;
        elseif count>1
            if ~isempty(linebuffer)
                if linebuffer(1)~=';'
                    remaining=linebuffer;
                    tok_cnt=1;
                    while ~isempty(remaining)
                        [str{tok_cnt},remaining] = strtok(remaining);
                        tok_cnt=tok_cnt+1;
                    end
                    for i=1:tok_cnt-1
                        if strcmp(str{i},'=')
                            ptr_name=i-1;
                            ptr_val=i+1;                    
                            network_file_section{count-1}.names=[network_file_section{count-1}.names; {str{ptr_name}}];
                            network_file_section{count-1}.values=[network_file_section{count-1}.values; {str{ptr_val}}];
                        end
                    end
                end
            end
        end
        line_number=line_number+1;
    end
    
    module_count=1;
    network_count=1;
    entry_count=1;
    cachegeometry_count=1;
    line_number=1;
    current=[];
    while ~feof(fid_memconfig)
        linebuffer=fgetl(fid_memconfig);
        if regexp(linebuffer,'\[Module .*\]')
            memconfig_module{module_count}.module_name=linebuffer(9:(length(linebuffer)-1));
            memconfig_module_name{module_count}=linebuffer(9:(length(linebuffer)-1));
            memconfig_module{module_count}.names=[];
            memconfig_module{module_count}.values=[];
            module_count=module_count+1;
            current=1;
        elseif regexp(linebuffer,'\[Network .*\]')
            memconfig_network{network_count}.network_name=linebuffer(10:(length(linebuffer)-1));
            memconfig_network_name{network_count}=linebuffer(10:(length(linebuffer)-1));
            memconfig_network{network_count}.names=[];
            memconfig_network{network_count}.values=[];
            network_count=network_count+1;
            current=2;
        elseif regexp(linebuffer,'\[Entry .*\]')
            memconfig_entry{entry_count}.entry_name=linebuffer(8:(length(linebuffer)-1));
            memconfig_entry_name{entry_count}=linebuffer(8:(length(linebuffer)-1));
            memconfig_entry{entry_count}.names=[];
            memconfig_entry{entry_count}.values=[];
            entry_count=entry_count+1;
            current=3;
        elseif regexp(linebuffer,'\[CacheGeometry .*\]')
            memconfig_cachegeometry{cachegeometry_count}.cachegeometry_name=linebuffer(16:(length(linebuffer)-1));
            memconfig_cachegeometry_name{cachegeometry_count}=linebuffer(16:(length(linebuffer)-1));
            memconfig_cachegeometry{cachegeometry_count}.names=[];
            memconfig_cachegeometry{cachegeometry_count}.values=[];
            cachegeometry_count=cachegeometry_count+1;
            current=4;
        elseif ~isempty(current)
            if ~isempty(linebuffer)
                if linebuffer(1)~=';'
                    remaining=linebuffer;
                    tok_cnt=1;
                    while ~isempty(remaining)
                        [str{tok_cnt},remaining] = strtok(remaining);
                        tok_cnt=tok_cnt+1;
                    end
                    for i=1:tok_cnt-1
                        if strcmp(str{i},'=')
                            ptr_name=i-1;
                            ptr_val=i+1;
                            if current==1
                                memconfig_module{module_count-1}.names=[memconfig_module{module_count-1}.names; {str{ptr_name}}];
                                memconfig_module{module_count-1}.values=[memconfig_module{module_count-1}.values; {str{ptr_val}}];
                            elseif current==2
                                memconfig_network{network_count-1}.names=[memconfig_network{network_count-1}.names; {str{ptr_name}}];
                                memconfig_network{network_count-1}.values=[memconfig_network{network_count-1}.values; {str{ptr_val}}];
                            elseif current==3
                                memconfig_entry{entry_count-1}.names=[memconfig_entry{entry_count-1}.names; {str{ptr_name}}];
                                memconfig_entry{entry_count-1}.values=[memconfig_entry{entry_count-1}.values; {str{ptr_val}}];                                
                            elseif current==4
                                memconfig_cachegeometry{cachegeometry_count-1}.names=[memconfig_cachegeometry{cachegeometry_count-1}.names; {str{ptr_name}}];
                                memconfig_cachegeometry{cachegeometry_count-1}.values=[memconfig_cachegeometry{cachegeometry_count-1}.values; {str{ptr_val}}];
                            else
                                display(['ERROR: while parsing the memconfig, the variable "current" took on the value ' num2str(current) ' but it should have a value in the set {[],1,2,3,4}'])
                            end
                        end
                    end
                end
            end
        end
        line_number=line_number+1;
    end   

    fclose(fid_mem);
    fclose(fid_pipeline);
    fclose(fid_memconfig);
    fclose(fid_network);

    %%%%%%%%%%%%%%%Generate all nessessary variables%%%%%%%%%%%%
    %num_cores
%     section_idx=find(strcmp(pipeline_file_section_name,'Config.General'));
%     property_idx=find(strcmp(pipeline_file_section{section_idx}.names,'Cores'));
%     num_cores=pipeline_file_section{section_idx}.values{property_idx}
    
    num_cores=str2num(reference(pipeline_file_section, pipeline_file_section_name, 'Config.General', 'Cores'));
    
    num_L2=num_cores;
    
    if num_L2>1
        homo_L2=0;
    else
        homo_L2=1;
    end
    
    if num_L3>0
        num_cache_level=3;
    elseif num_L2>0
        num_cache_level=2;
    else
        num_cache_level=1;
    end
    
    %check if num_cores aggress with length(memconfig_entry_name)
    if num_cores~=length(memconfig_entry_name)
        display(['WARNING: num_cores=' num2str(num_cores) ' but length(memconfig_entry_name)=' num2str(length(memconfig_entry_name))])
    end
    
    cycles=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Global','Cycles'));
    num_threads=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.General','Threads'));
    decode_width=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Pipeline','DecodeWidth'));
    issue_width=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Pipeline','IssueWidth'));
    commit_width=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Pipeline','CommitWidth'));
    alu_per_core=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.FunctionalUnits','IntAdd.Count'));
    mul_per_core=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.FunctionalUnits','IntMult.Count'));
    fpu_per_core=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.FunctionalUnits','FloatSimple.Count'));
    size_IQ=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','IqSize'));
    size_ROB=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','RobSize'));
    size_IRF=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','RfIntSize'));
    size_FRF=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','RfFpSize'));
    size_LSQ=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','LsqSize'));
    size_RAS=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','RAS.Size'));
    
    size_fetch_queue=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.Queues','FetchQueueSize'));
    size_bpred=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','TwoLevel.L1Size'));
    size_bimod=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','Bimod.Size'));
    size_hist=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','TwoLevel.HistorySize'));
    
    for core_number=1:num_cores
        dcache_ports(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'Ports'));
        icache_ports(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[icache_mod num2str(core_number-1)],'Ports'));
        L2_ports(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Ports'));
        total_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Total'));
        int_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Integer'));
        fp_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.FloatingPoint'));
        branch_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Ctrl'));
        branch_mispred(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Commit.Mispred'));
        load_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Uop.load'));
        store_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Uop.store'));
        committed_inst(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Commit.Total'));
        committed_int(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Commit.Integer'));
        committed_fp(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Commit.FloatingPoint'));
        pipeline_duty(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Commit.DutyCycle'));
        ROB_reads(core_number)=0;
        ROB_writes(core_number)=0;
        RAT_reads(core_number)=0;
        RAT_writes(core_number)=0;
        RAT_fp_reads(core_number)=0;
        RAT_fp_writes(core_number)=0;
        IQ_reads(core_number)=0;
        IQ_writes(core_number)=0;
        IQ_wakeup(core_number)=0;
        ireg_reads(core_number)=0;
        freg_reads(core_number)=0;
        ireg_writes(core_number)=0;
        freg_writes(core_number)=0;
        for thread_number=1:num_threads
            ROB_reads(core_number)=ROB_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'ROB.Reads'));
            ROB_writes(core_number)=ROB_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'ROB.Writes'));
            RAT_reads(core_number)=RAT_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.IntReads'))+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.FpReads'));
            RAT_writes(core_number)=RAT_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.IntWrites'))+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.FpWrites'));
            RAT_fp_reads(core_number)=RAT_fp_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.FpReads'));
            RAT_fp_writes(core_number)=RAT_fp_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RAT.FpWrites'));
            IQ_reads(core_number)=IQ_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'IQ.Reads'));
            IQ_writes(core_number)=IQ_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'IQ.Writes'));
            %IQ_wakeup(core_number)=IQ_wakeup(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'IQ.WakeupAccesses'));
            ireg_reads(core_number)=ireg_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RF_Int.Reads'));
            freg_reads(core_number)=freg_reads(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RF_Fp.Reads'));
            ireg_writes(core_number)=ireg_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RF_Int.Writes'));
            freg_writes(core_number)=freg_writes(core_number)+str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1) ' Thread ' num2str(thread_number-1)],'RF_Fp.Writes'));
        end
        func_calls(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.Uop.call'));
        ctx_switches(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Dispatch.WndSwitch'));
        ialu_accesses(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Issue.Integer'));
        fpu_accesses(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Issue.FloatingPoint'));
        mul_accesses(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,['Core ' num2str(core_number-1)],'Issue.Uop.mult'));
        icache_capacity(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'Sets'))*str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'Assoc'))*str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'BlockSize'));
        icache_blockwidth(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'BlockSize'));
        icache_assoc(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'Assoc'));
        icache_latency(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,icache_geo,'Latency'));
        icache_read_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[icache_mod num2str(core_number-1)],'BlockingReads'));
        icache_read_misses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[icache_mod num2str(core_number-1)],'ReadMisses'));
        icache_conflicts(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[icache_mod num2str(core_number-1)],'Evictions'));
        dcache_capacity(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'Sets'))*str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'Assoc'))*str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'BlockSize'));
        dcache_blockwidth(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'BlockSize'));
        dcache_assoc(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'Assoc'));
        dcache_latency(core_number)=str2num(reference(memconfig_cachegeometry,memconfig_cachegeometry_name,dcache_geo,'Latency'));
        dcache_read_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'Reads'));
%         dcache_read_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'BlockingReads'));
        dcache_read_misses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'ReadMisses'));
        dcache_write_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'Writes'));
%         dcache_write_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'BlockingWrites'));
        dcache_write_misses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'WriteMisses'));
        dcache_conflicts(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[dcache_mod num2str(core_number-1)],'Evictions'));
        BTB_capacity(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','BTB.Sets'))*str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','BTB.Assoc'))*4;
        BTB_assoc(core_number)=str2num(reference(pipeline_file_section,pipeline_file_section_name,'Config.BranchPredictor','BTB.Assoc'));
    end
    total_mem_accesses=0;
    total_mem_reads=0;
    total_mem_writes=0;
    for i=1:num_MC
        mem_latency{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'DataLatency'));
        mem_blocksize{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'BlockSize'));
        mem_accesses{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'Accesses'));
        mem_reads{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'Reads'));
        mem_writes{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'Writes'));
        mem_ports{i}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{i},'Ports'));
        total_mem_accesses=total_mem_accesses+mem_accesses{i};
        total_mem_reads=total_mem_reads+mem_reads{i};
        total_mem_writes=total_mem_writes+mem_writes{i};
    end
%     mem_latency{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{2},'Latency'));
%     mem_blocksize{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{1},'BlockSize'));
%     mem_accesses{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{2},'Accesses'));
%     mem_reads{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{2},'Reads'));
%     mem_writes{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{2},'Writes'));
%     mem_ports{2}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{2},'Ports'));
%     mem_latency{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'Latency'));
%     mem_blocksize{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'BlockSize'));
%     mem_accesses{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'Accesses'));
%     mem_reads{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'Reads'));
%     mem_writes{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'Writes'));
%     mem_ports{3}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{3},'Ports'));
%     mem_latency{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'Latency'));
%     mem_blocksize{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'BlockSize'));
%     mem_accesses{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'Accesses'));
%     mem_reads{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'Reads'));
%     mem_writes{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'Writes'));
%     mem_ports{4}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{4},'Ports'));
%     mem_latency{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'Latency'));
%     mem_blocksize{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'BlockSize'));
%     mem_accesses{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'Accesses'));
%     mem_reads{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'Reads'));
%     mem_writes{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'Writes'));
%     mem_ports{5}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{5},'Ports'));
%     mem_latency{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'Latency'));
%     mem_blocksize{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'BlockSize'));
%     mem_accesses{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'Accesses'));
%     mem_reads{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'Reads'));
%     mem_writes{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'Writes'));
%     mem_ports{6}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{6},'Ports'));
%     mem_latency{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'Latency'));
%     mem_blocksize{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'BlockSize'));
%     mem_accesses{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'Accesses'));
%     mem_reads{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'Reads'));
%     mem_writes{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'Writes'));
%     mem_ports{7}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{7},'Ports'));
%     mem_latency{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'Latency'));
%     mem_blocksize{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'BlockSize'));
%     mem_accesses{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'Accesses'));
%     mem_reads{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'Reads'));
%     mem_writes{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'Writes'));
%     mem_ports{8}=str2num(reference(mem_file_section,mem_file_section_name,mem_mod{8},'Ports'));
    %%%%%%%%%%FIXME%%%%%%%%%
    network_accesses=str2num(reference(network_file_section,network_file_section_name,'Network.net-l1-l2','Transfers'));
    network_duty=network_accesses/cycles;
    %%%%%%%%%%%%%%%%%%%%%%%% Multiple L2s
    for core_number=1:num_cores
        L2_capacity(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Sets'))*str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Ways'))*str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'BlockSize'));
        L2_blockwidth(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'BlockSize'));
        L2_assoc(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Ways'));
        L2_latency(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'DataLatency'));
        L2_read_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Reads'));
        L2_write_accesses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Writes'));
        L2_read_misses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'ReadMisses'));
        L2_write_misses(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'WriteMisses'));
        L2_conflicts(core_number)=str2num(reference(mem_file_section,mem_file_section_name,[L2_mod num2str(core_number-1)],'Evictions'));
        L2_duty(core_number)=1;
    end
    
%     %count number of L2s
%     L2_modules=[];
%     MM_modules=[];
%     for i=1:length(memconfig_entry_name)
%         property_idx=find(strcmp(memconfig_entry{i}.names,'DataModule'));
%         datamodule{i}=memconfig_entry{i}.values{property_idx};
%         property_idx=find(strcmp(memconfig_entry{i}.names,'InstModule'));
%         instmodule{i}=memconfig_entry{i}.values{property_idx};
%         
%         section_idx=find(strcmp(memconfig_module_name,datamodule{i}));
%         property_idx=find(strcmp(memconfig_module{section_idx}.names,'LowModules'));
%         
%         LowModules=reference(memconfig_module,memconfig_module_name,datamodule{i},'LowModules');
%         
%         if isempty(find(strcmp(L2_modules,LowModules)))
% %         if isempty(find(strcmp(L2_modules,memconfig_module{section_idx}.values{property_idx})))
%             %check if this module is actualy MM or L2
% %             sec_id=find(strcmp(memconfig_module_name,memconfig_module{section_idx}.values{property_idx}));
% %             prop_id=find(strcmp(memconfig_module{sec_id}.names,'Type'));
% %             if strcmp(memconfig_module{sec_id}.values{prop_id},'Cache')            
% %                 L2_modules=[L2_modules; {memconfig_module{section_idx}.values{property_idx}}];
% %             else
% %                 MM_modules=[MM_modules; {memconfig_module{section_idx}.values{property_idx}}];
% %             end
%             
%             type=reference(memconfig_module, memconfig_module_name, memconfig_module{section_idx}.values{property_idx}, 'Type');
%             if strcmp(type,'Cache')
%                 L2_modules=[L2_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             else
%                 MM_modules=[MM_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             end
%         end
%         
%         section_idx=find(strcmp(memconfig_module_name,instmodule{i}));
%         property_idx=find(strcmp(memconfig_module{section_idx}.names,'LowModules'));
%         if isempty(find(strcmp(L2_modules,memconfig_module{section_idx}.values{property_idx})))
%             %check if this module is actualy MM or L2
%             sec_id=find(strcmp(memconfig_module_name,memconfig_module{section_idx}.values{property_idx}));
%             prop_id=find(strcmp(memconfig_module{sec_id}.names,'Type'));
%             if strcmp(memconfig_module{sec_id}.values{prop_id},'Cache')
%                 L2_modules=[L2_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             else
%                 MM_modules=[MM_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             end
%         end
%     end
%     L2_modules
%     MM_modules
%     num_L2=length(L2_modules);
%     num_MM=length(MM_modules);
%     
%     %count number of L3s
%     L3_modules=[];
%     for i=1:num_L2
%         section_idx=find(strcmp(memconfig_module_name,L2_modules{i}));
%         property_idx=find(strcmp(memconfig_module{section_idx}.names,'LowModules'));
%         if isempty(find(strcmp(L3_modules,memconfig_module{section_idx}.values{property_idx})))
%             %check if this module is actualy MM or L3
%             sec_id=find(strcmp(memconfig_module_name,memconfig_module{section_idx}.values{property_idx}));
%             prop_id=find(strcmp(memconfig_module{sec_id}.names,'Type'));
%             if strcmp(memconfig_module{sec_id}.values{prop_id},'Cache')
%                 L3_modules=[L3_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             else
%                 MM_modules=[MM_modules; {memconfig_module{section_idx}.values{property_idx}}];
%             end
%         end
%     end
%     L3_modules
%     MM_modules
%     num_L3=length(L3_modules);
%     num_MM=length(MM_modules);
%     
%     
%     %Make sure there is not L4

    %
    
    %%%%%%%%%%%%%%%%%%%%%%%Start generateing XML%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fid_xml=fopen(xml_file,'w');
    
    fprintf(fid_xml,'<?xml version="1.0" ?>\n');
    fprintf(fid_xml,'<component id="root" name="root">\n');
    fprintf(fid_xml,'	<component id="system" name="system">\n');
    fprintf(fid_xml,'		<!--McPAT will skip the components if number is set to 0 -->\n');
    fprintf(fid_xml,['		<param name="number_of_cores" value="' num2str(num_cores) '"/>\n']);
    if L1_dir==1
        fprintf(fid_xml,['		<param name="number_of_L1Directories" value="' num2str(num_cores) '"/>\n']);
    else
        fprintf(fid_xml,['		<param name="number_of_L1Directories" value="0"/>\n']);
    end
    if L2_dir==1
        fprintf(fid_xml,['		<param name="number_of_L2Directories" value="' num2str(num_cores) '"/>\n']);
    else
        fprintf(fid_xml,['		<param name="number_of_L2Directories" value="0"/>\n']);
    end
    fprintf(fid_xml,['		<param name="number_of_L2s" value="' num2str(num_L2) '"/> <!-- This number means how many L2 clusters in each cluster there can be multiple banks/ports -->\n']);
    fprintf(fid_xml,['		<param name="number_of_L3s" value="' num2str(num_L3) '"/> <!-- This number means how many L3 clusters -->\n']);
    fprintf(fid_xml,['		<param name="number_of_NoCs" value="1"/>\n']);
    fprintf(fid_xml,['		<param name="homogeneous_cores" value="' num2str(homo_core) '"/><!--1 means homo -->\n']);
    fprintf(fid_xml,['		<param name="homogeneous_L2s" value="' num2str(homo_L2) '"/>\n']);
    fprintf(fid_xml,['		<param name="homogeneous_L1Directorys" value="' num2str(homo_L1_dir) '"/>\n']);
    fprintf(fid_xml,['		<param name="homogeneous_L2Directorys" value="' num2str(homo_L2_dir) '"/>\n']);
    fprintf(fid_xml,['		<param name="homogeneous_L3s" value="' num2str(homo_L3) '"/>\n']);
    fprintf(fid_xml,['		<param name="homogeneous_ccs" value="1"/><!--cache coherece hardware -->\n']);
    fprintf(fid_xml,['		<param name="homogeneous_NoCs" value="1"/>\n']);
    fprintf(fid_xml,['		<param name="core_tech_node" value="' num2str(tech) '"/><!-- nm -->\n']);
    fprintf(fid_xml,['		<param name="target_core_clockrate" value="' num2str(clk) '"/><!--MHz -->\n']);
    fprintf(fid_xml,['		<param name="temperature" value="' num2str(temp) '"/> <!-- Kelvin -->\n']);
    fprintf(fid_xml,['		<param name="number_cache_levels" value="' num2str(num_cache_level) '"/>\n']);
    fprintf(fid_xml,['		<param name="interconnect_projection_type" value="0"/><!--0: agressive wire technology; 1: conservative wire technology -->\n']);
    fprintf(fid_xml,['		<param name="device_type" value="0"/><!--0: HP(High Performance Type); 1: LSTP(Low standby power) 2: LOP (Low Operating Power)  -->\n']);
    fprintf(fid_xml,['		<param name="longer_channel_device" value="1"/><!-- 0 no use; 1 use when possible -->\n']);
    fprintf(fid_xml,['		<param name="machine_bits" value="32"/>\n']);
    fprintf(fid_xml,['		<param name="virtual_address_width" value="32"/>\n']);
    fprintf(fid_xml,['		<param name="physical_address_width" value="32"/>\n']);
    fprintf(fid_xml,['    		<param name="virtual_memory_page_size" value="4096"/>\n']);
    fprintf(fid_xml,['		<stat name="total_cycles" value="' num2str(cycles) '"/>\n']);
    fprintf(fid_xml,['		<stat name="idle_cycles" value="0"/>\n']);
    fprintf(fid_xml,['		<stat name="busy_cycles"  value="' num2str(cycles) '"/>\n']);
    fprintf(fid_xml,['			<!--This page size(B) is complete different from the page size in Main memo secction. this page size is the size of \n']);
    fprintf(fid_xml,['			virtual memory from OS/Archi perspective; the page size in Main memo secction is the actuall physical line in a DRAM bank  -->\n']);
    fprintf(fid_xml,['		<!-- *********************** cores ******************* -->\n']);
    for core_number=1:num_cores
        fprintf(fid_xml,['		<component id="system.core' num2str(core_number-1) '" name="core' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['			<!-- Core property -->\n']);
        fprintf(fid_xml,['			<param name="clock_rate" value="' num2str(clk) '"/>\n']);
        fprintf(fid_xml,['			<param name="instruction_length" value="32"/>\n']);
        fprintf(fid_xml,['			<param name="opcode_width" value="9"/>\n']);
        fprintf(fid_xml,['			<!-- address width determins the tag_width in Cache, LSQ and buffers in cache controller \n']);
        fprintf(fid_xml,['			default value is machine_bits, if not set --> \n']);
        fprintf(fid_xml,['			<param name="machine_type" value="0"/><!-- 1 inorder; 0 OOO-->\n']);
        fprintf(fid_xml,['			<!-- inorder/OoO -->\n']);
        fprintf(fid_xml,['			<param name="number_hardware_threads" value="' num2str(num_threads) '"/>\n']);
        fprintf(fid_xml,['			<!-- number_instruction_fetch_ports(icache ports) is always 1 in single-thread processor,\n']);
        fprintf(fid_xml,['			it only may be more than one in SMT processors. BTB ports always equals to fetch ports since \n']);
        fprintf(fid_xml,['			branch information in consective branch instructions in the same fetch group can be read out from BTB once.--> \n']);
        fprintf(fid_xml,['			<param name="fetch_width" value="' num2str(decode_width) '"/>\n']);
        fprintf(fid_xml,['			<!-- fetch_width determins the size of cachelines of L1 cache block -->\n']);
        fprintf(fid_xml,['			<param name="number_instruction_fetch_ports" value="' num2str(decode_width) '"/>\n']);
        fprintf(fid_xml,['			<param name="decode_width" value="' num2str(decode_width) '"/>\n']);
        fprintf(fid_xml,['			<!-- decode_width determins the number of ports of the \n']);
        fprintf(fid_xml,['			renaming table (both RAM and CAM) scheme -->\n']);
        fprintf(fid_xml,['			<param name="issue_width" value="' num2str(issue_width) '"/>\n']);
        fprintf(fid_xml,['			<!-- issue_width determins the number of ports of Issue window and other logic \n']);
        fprintf(fid_xml,['			as in the complexity effective proccessors paper; issue_width==dispatch_width -->\n']);
        fprintf(fid_xml,['			<param name="commit_width" value="' num2str(commit_width) '"/>\n']);
        fprintf(fid_xml,['			<!-- commit_width determins the number of ports of register files -->\n']);
        fprintf(fid_xml,['			<param name="fp_issue_width" value="' num2str(issue_width) '"/>\n']);
        fprintf(fid_xml,['			<param name="prediction_width" value="' num2str(decode_width) '"/> \n']);
        fprintf(fid_xml,['			<!-- number of branch instructions can be predicted simultannouesl-->\n']);
        fprintf(fid_xml,['			<!-- Current version of McPAT does not distinguish int and floating point pipelines \n']);
        fprintf(fid_xml,['			Theses parameters are reserved for future use.--> \n']);
        fprintf(fid_xml,['			<param name="pipelines_per_core" value="1,0"/>\n']);
        fprintf(fid_xml,['			<!--integer_pipeline and floating_pipelines, if the floating_pipelines is 0, then the pipeline is shared-->\n']);
        fprintf(fid_xml,['			<param name="pipeline_depth" value="6,6"/>\n']);
        fprintf(fid_xml,['			<!-- pipeline depth of int and fp, if pipeline is shared, the second number is the average cycles of fp ops -->\n']);
        fprintf(fid_xml,['			<!-- issue and exe unit-->\n']);
        fprintf(fid_xml,['			<param name="ALU_per_core" value="' num2str(alu_per_core) '"/>\n']);
        fprintf(fid_xml,['			<!-- contains an adder, a shifter, and a logical unit -->\n']);
        fprintf(fid_xml,['			<param name="MUL_per_core" value="' num2str(mul_per_core) '"/>\n']);
        fprintf(fid_xml,['			<!-- For MUL and Div -->\n']);
        fprintf(fid_xml,['			<param name="FPU_per_core" value="' num2str(fpu_per_core) '"/>	\n']);
        fprintf(fid_xml,['			<!-- buffer between IF and ID stage -->\n']);
        fprintf(fid_xml,['			<param name="instruction_buffer_size" value="' num2str(size_fetch_queue) '"/>\n']);
        fprintf(fid_xml,['			<!-- buffer between ID and sche/exe stage -->\n']);
%         fprintf(fid_xml,['			<param name="decoded_stream_buffer_size" value="0"/>\n']);
        fprintf(fid_xml,['			<param name="instruction_window_scheme" value="0"/><!-- 0 PHYREG based, 1 RSBASED-->\n']);
        fprintf(fid_xml,['			<!-- McPAT support 2 types of OoO cores, RS based and physical reg based-->\n']);
        fprintf(fid_xml,['			<param name="instruction_window_size" value="' num2str(size_IQ/2) '"/>\n']);
        fprintf(fid_xml,['			<param name="fp_instruction_window_size" value="' num2str(size_IQ/2) '"/>\n']);
        fprintf(fid_xml,['			<!-- the instruction issue Q as in Alpha 21264; The RS as in Intel P6 -->\n']);
        fprintf(fid_xml,['			<param name="ROB_size" value="' num2str(size_ROB) '"/>\n']);
        fprintf(fid_xml,['			<!-- each in-flight instruction has an entry in ROB -->\n']);
        fprintf(fid_xml,['			<!-- registers -->\n']);
        fprintf(fid_xml,['			<param name="archi_Regs_IRF_size" value="' num2str(size_arch_IRF) '"/>\n']);
        fprintf(fid_xml,['			<param name="archi_Regs_FRF_size" value="' num2str(size_arch_FRF) '"/>\n']);
        fprintf(fid_xml,['			<!--  if OoO processor, phy_reg number is needed for renaming logic, \n']);
        fprintf(fid_xml,['			renaming logic is for both integer and floating point insts.  -->\n']);
        fprintf(fid_xml,['			<param name="phy_Regs_IRF_size" value="' num2str(size_IRF) '"/>\n']);
        fprintf(fid_xml,['			<param name="phy_Regs_FRF_size" value="' num2str(size_FRF) '"/>\n']);
        fprintf(fid_xml,['			<!-- rename logic -->\n']);
        fprintf(fid_xml,['			<param name="rename_scheme" value="0"/>\n']);
        fprintf(fid_xml,['			<!-- can be RAM based(0) or CAM based(1) rename scheme \n']);
        fprintf(fid_xml,['			RAM-based scheme will have free list, status table;\n']);
        fprintf(fid_xml,['			CAM-based scheme have the valid bit in the data field of the CAM \n']);
        fprintf(fid_xml,['			both RAM and CAM need RAM-based checkpoint table, checkpoint_depth=# of in_flight instructions;\n']);
        fprintf(fid_xml,['			Detailed RAT Implementation see TR -->\n']);
%         fprintf(fid_xml,['			<param name="register_windows_size" value="8"/>\n']);
        fprintf(fid_xml,['			<!-- how many windows in the windowed register file, sun processors;\n']);
        fprintf(fid_xml,['			no register windowing is used when this number is 0 -->\n']);
        fprintf(fid_xml,['			<!-- In OoO cores, loads and stores can be issued whether inorder(Pentium Pro) or (OoO)out-of-order(Alpha),\n']);
        fprintf(fid_xml,['			They will always try to exeute out-of-order though. -->\n']);
        fprintf(fid_xml,['			<param name="LSU_order" value="OoO"/>\n']);
        fprintf(fid_xml,['			<param name="store_buffer_size" value="' num2str(size_LSQ/2) '"/>\n']);
        fprintf(fid_xml,['			<!-- By default, in-order cores do not have load buffers -->\n']);
        fprintf(fid_xml,['			<param name="load_buffer_size" value="' num2str(size_LSQ/2) '"/>	\n']);
        fprintf(fid_xml,['			<!-- number of ports refer to sustainable concurrent memory accesses --> \n']);
        fprintf(fid_xml,['			<param name="memory_ports" value="' num2str(dcache_ports(core_number)) '"/>	\n']);
%         fprintf(fid_xml,['			<param name="memory_ports" value="' num2str(n*mem_ports{1}) '"/>	\n']);
        fprintf(fid_xml,['			<!-- max_allowed_in_flight_memo_instructions determins the # of ports of load and store buffer\n']);
        fprintf(fid_xml,['			as well as the ports of Dcache which is connected to LSU -->	\n']);
        fprintf(fid_xml,['			<!-- dual-pumped Dcache can be used to save the extra read/write ports -->\n']);
        fprintf(fid_xml,['			<param name="RAS_size" value="' num2str(size_RAS) '"/>\n']);
        fprintf(fid_xml,['			<!-- general stats, defines simulation periods;require total, idle, and busy cycles for senity check  -->\n']);
        fprintf(fid_xml,['			<!-- please note: if target architecture is X86, then all the instrucions refer to (fused) micro-ops -->\n']);
        fprintf(fid_xml,['			<stat name="total_instructions" value="' num2str(total_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="int_instructions" value="' num2str(int_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="fp_instructions" value="' num2str(fp_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="branch_instructions" value="' num2str(branch_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="branch_mispredictions" value="' num2str(branch_mispred(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="load_instructions" value="' num2str(load_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="store_instructions" value="' num2str(store_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="committed_instructions" value="' num2str(committed_inst(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="committed_int_instructions" value="' num2str(committed_int(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="committed_fp_instructions" value="' num2str(committed_fp(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="pipeline_duty_cycle" value="' num2str(pipeline_duty(core_number)) '"/><!--<=1, runtime_ipc/peak_ipc; averaged for all cores if homogenous -->\n']);
        fprintf(fid_xml,['			<!-- the following cycle stats are used for heterogeneouse cores only, \n']);  
        fprintf(fid_xml,['				please ignore them if homogeneouse cores -->\n']);
        fprintf(fid_xml,['			<stat name="total_cycles" value="' num2str(cycles) '"/>\n']);
        fprintf(fid_xml,['		    <stat name="idle_cycles" value="0"/>\n']);
        fprintf(fid_xml,['		    <stat name="busy_cycles"  value="' num2str(cycles) '"/>\n']);
        fprintf(fid_xml,['			<!-- instruction buffer stats -->\n']);
        fprintf(fid_xml,['			<!-- ROB stats, both RS and Phy based OoOs have ROB\n']);
        fprintf(fid_xml,['			performance simulator should capture the difference on accesses,\n']);
        fprintf(fid_xml,['			otherwise, McPAT has to guess based on number of commited instructions. -->\n']);
        fprintf(fid_xml,['			<stat name="ROB_reads" value="' num2str(ROB_reads(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="ROB_writes" value="' num2str(ROB_writes(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<!-- RAT accesses -->\n']);
        fprintf(fid_xml,['			<stat name="rename_reads" value="' num2str(RAT_reads(core_number)) '"/>\n']);  
        fprintf(fid_xml,['			<stat name="rename_writes" value="' num2str(RAT_writes(core_number)) '"/>\n']);  
        fprintf(fid_xml,['			<stat name="fp_rename_reads" value="' num2str(RAT_fp_reads(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="fp_rename_writes" value="' num2str(RAT_fp_writes(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<!-- decode and rename stage use this, should be total ic - nop -->\n']);
        fprintf(fid_xml,['			<!-- Inst window stats -->\n']);
        fprintf(fid_xml,['			<stat name="inst_window_reads" value="' num2str(IQ_reads(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="inst_window_writes" value="' num2str(IQ_writes(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="inst_window_wakeup_accesses" value="' num2str(IQ_wakeup(core_number)) '"/>\n']);
%         fprintf(fid_xml,['			<stat name="fp_inst_window_reads" value="' num2str(IQ_fp_reads(core_number)) '"/>\n']);
%         fprintf(fid_xml,['			<stat name="fp_inst_window_writes" value="' num2str(IQ_fp_writes(core_number)) '"/>\n']);
%         fprintf(fid_xml,['			<stat name="fp_inst_window_wakeup_accesses" value="' num2str(IQ_fp_wakeup(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<!--  RF accesses -->\n']);
        fprintf(fid_xml,['			<stat name="int_regfile_reads" value="' num2str(ireg_reads(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="float_regfile_reads" value="' num2str(freg_reads(core_number)) '"/>\n']);  
        fprintf(fid_xml,['			<stat name="int_regfile_writes" value="' num2str(ireg_writes(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="float_regfile_writes" value="' num2str(freg_writes(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<!-- accesses to the working reg -->\n']);
        fprintf(fid_xml,['			<stat name="function_calls" value="' num2str(func_calls(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="context_switches" value="' num2str(ctx_switches(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<!-- Number of Windowes switches (number of function calls and returns)-->\n']);
        fprintf(fid_xml,['			<!-- Alu stats by default, the processor has one FPU that includes the divider and \n']);
        fprintf(fid_xml,['			 multiplier. The fpu accesses should include accesses to multiplier and divider  -->\n']);
        fprintf(fid_xml,['			<stat name="ialu_accesses" value="' num2str(ialu_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="fpu_accesses" value="' num2str(fpu_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['			<stat name="mul_accesses" value="' num2str(mul_accesses(core_number)) '"/>\n']);
%         fprintf(fid_xml,['			<stat name="cdb_alu_accesses" value="' num2str() '"/>\n']);  
%         fprintf(fid_xml,['			<stat name="cdb_mul_accesses" value="' num2str() '"/>\n']);
%         fprintf(fid_xml,['			<stat name="cdb_fpu_accesses" value="' num2str() '"/>\n']);
        fprintf(fid_xml,['			<!-- multiple cycle accesses should be counted multiple times, \n']);
        fprintf(fid_xml,['			otherwise, McPAT can use internal counter for different floating point instructions \n']);
        fprintf(fid_xml,['			to get final accesses. But that needs detailed info for floating point inst mix -->\n']);
        fprintf(fid_xml,['			<!--  currently the performance simulator should \n']);
        fprintf(fid_xml,['			make sure all the numbers are final numbers, \n']);
        fprintf(fid_xml,['			including the explicit read/write accesses, \n']);
        fprintf(fid_xml,['			and the implicite accesses such as replacements and etc.\n']);
        fprintf(fid_xml,['			Future versions of McPAT may be able to reason the implicite access\n']);
        fprintf(fid_xml,['			based on param and stats of last level cache\n']);
        fprintf(fid_xml,['			The same rule applies to all cache access stats too!  -->\n']);  
%         fprintf(fid_xml,['			<!-- following is AF for max power computation. \n']);
%         fprintf(fid_xml,['				Do not change them, unless you understand them-->\n']);
%         fprintf(fid_xml,['			<stat name="IFU_duty_cycle" value="0.5"/>\n']);
%         fprintf(fid_xml,['			<stat name="LSU_duty_cycle" value="0.25"/>\n']);
%         fprintf(fid_xml,['			<stat name="MemManU_I_duty_cycle" value="0.5"/>\n']);
%         fprintf(fid_xml,['			<stat name="MemManU_D_duty_cycle" value="0.25"/>\n']);
%         fprintf(fid_xml,['			<stat name="ALU_duty_cycle" value="0.9"/>\n']);
%         fprintf(fid_xml,['			<stat name="MUL_duty_cycle" value="0.75"/>\n']);
%         fprintf(fid_xml,['			<stat name="FPU_duty_cycle" value="0.6"/>\n']);
% %         fprintf(fid_xml,['			<!--FPU also handles Mul/div -->\n']);
%         fprintf(fid_xml,['			<stat name="ALU_cdb_duty_cycle" value="0.9"/>\n']);
%         fprintf(fid_xml,['			<stat name="MUL_cdb_duty_cycle" value="0.75"/>\n']);  
%         fprintf(fid_xml,['			<stat name="FPU_cdb_duty_cycle" value="0.6"/>	\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.predictor" name="PBT' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<!-- branch predictor; tournament predictor see Alpha implementation -->\n']);
        fprintf(fid_xml,['				<param name="local_predictor_size" value="' num2str(size_hist) ',1"/>\n']);
        fprintf(fid_xml,['				<param name="local_predictor_entries" value="' num2str(size_bpred) '"/>\n']);
        fprintf(fid_xml,['				<param name="global_predictor_entries" value="' num2str(size_bimod) '"/>\n']);
        fprintf(fid_xml,['				<param name="global_predictor_bits" value="2"/>\n']);
        fprintf(fid_xml,['				<param name="chooser_predictor_entries" value="' num2str(size_bimod) '"/>\n']);
        fprintf(fid_xml,['				<param name="chooser_predictor_bits" value="2"/>\n']);
        fprintf(fid_xml,['				<!-- These parameters can be combined like below in next version\n']);  
        fprintf(fid_xml,['				<param name="load_predictor" value="10,3,1024"/>\n']);
        fprintf(fid_xml,['				<param name="global_predictor" value="4096,2"/>\n']);
        fprintf(fid_xml,['				<param name="predictor_chooser" value="4096,2"/>\n']);
        fprintf(fid_xml,['				-->\n']);
        fprintf(fid_xml,['			</component>\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.itlb" name="itlb' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<param name="number_entries" value="64"/>\n']);
        fprintf(fid_xml,['				<stat name="total_accesses" value="0"/>\n']);
        fprintf(fid_xml,['				<stat name="total_misses" value="0"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="0"/>	\n']);
        fprintf(fid_xml,['				<!-- there is no write requests to itlb although writes happen to itlb after miss, \n']);
        fprintf(fid_xml,['				which is actually a replacement -->\n']);  
        fprintf(fid_xml,['			</component>\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.icache" name="icache' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<!-- there is no write requests to itlb although writes happen to it after miss, \n']);
        fprintf(fid_xml,['				which is actually a replacement -->\n']);
        fprintf(fid_xml,['				<param name="icache_config" value="' num2str(icache_capacity(core_number)) ',' num2str(icache_blockwidth(core_number)) ',' num2str(icache_assoc(core_number)) ',1,' num2str(icache_latency(core_number)) ',' num2str(icache_latency(core_number)) ',8,1"/>\n']);
        fprintf(fid_xml,['				<!-- the parameters are capacity,block_width, associativity, bank, throughput w.r.t. core clock, latency w.r.t. core clock,output_width, cache policy -->\n']);
        fprintf(fid_xml,['				<!-- cache_policy;//0 no write or write-though with non-write allocate;1 write-back with write-allocate -->\n']);
        fprintf(fid_xml,['				<param name="buffer_sizes" value="' num2str(icache_blockwidth(core_number)) ', ' num2str(icache_blockwidth(core_number)) ', ' num2str(icache_blockwidth(core_number)) ',0"/>\n']);
        fprintf(fid_xml,['				<!-- cache controller buffer sizes: miss_buffer_size(MSHR),fill_buffer_size,prefetch_buffer_size,wb_buffer_size--> \n']);
        fprintf(fid_xml,['				<stat name="read_accesses" value="' num2str(icache_read_accesses(core_number)) '"/>\n']);  
        fprintf(fid_xml,['				<stat name="read_misses" value="' num2str(icache_read_misses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="' num2str(icache_conflicts(core_number)) '"/>	\n']);
        fprintf(fid_xml,['			</component>\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.dtlb" name="dtlb' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<param name="number_entries" value="128"/>\n']);
        fprintf(fid_xml,['				<stat name="total_accesses" value="0"/>\n']);
        fprintf(fid_xml,['				<stat name="total_misses" value="0"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="0"/>	\n']);
        fprintf(fid_xml,['			</component>\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.dcache" name="dcache' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['			        <!-- all the buffer related are optional -->\n']);
        fprintf(fid_xml,['				<param name="dcache_config" value="' num2str(dcache_capacity(core_number)) ',' num2str(dcache_blockwidth(core_number)) ',' num2str(dcache_assoc(core_number)) ',1, ' num2str(dcache_latency(core_number)) ',' num2str(dcache_latency(core_number)) ', 8,1"/>\n']);  
        fprintf(fid_xml,['				<param name="buffer_sizes" value="' num2str(dcache_blockwidth(core_number)) ', ' num2str(dcache_blockwidth(core_number)) ', ' num2str(dcache_blockwidth(core_number)) ', ' num2str(dcache_blockwidth(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<!-- cache controller buffer sizes: miss_buffer_size(MSHR),fill_buffer_size,prefetch_buffer_size,wb_buffer_size-->	\n']);
        fprintf(fid_xml,['				<stat name="read_accesses" value="' num2str(dcache_read_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="write_accesses" value="' num2str(dcache_write_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="read_misses" value="' num2str(dcache_read_misses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="write_misses" value="' num2str(dcache_write_misses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="' num2str(dcache_conflicts(core_number)) '"/>	\n']);
        fprintf(fid_xml,['			</component>\n']);
        fprintf(fid_xml,['			<component id="system.core' num2str(core_number-1) '.BTB" name="BTB' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['			        <!-- all the buffer related are optional -->\n']);  
        fprintf(fid_xml,['				<param name="BTB_config" value="' num2str(BTB_capacity(core_number)) ',4,' num2str(BTB_assoc(core_number)) ',1, 1,1"/>\n']);
        fprintf(fid_xml,['				<!-- the parameters are capacity,block_width,associativity,bank, throughput w.r.t. core clock, latency w.r.t. core clock,-->\n']);
        fprintf(fid_xml,['			</component>\n']); 
        fprintf(fid_xml,['	    </component>\n']); %end of core component
    end
    for core_number=1:1
        fprintf(fid_xml,['		<component id="system.L1Directory' num2str(core_number-1) '" name="L1Directory' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<param name="Directory_type" value="1"/>\n']);
        fprintf(fid_xml,['			    <!--0 cam based shadowed tag. 1 directory cache -->\n']);
        fprintf(fid_xml,['				<param name="Dir_config" value="1024,2,0,1,1,1, 8"/>\n']);
        fprintf(fid_xml,['				<!-- the parameters are capacity,block_width, associativity,bank, throughput w.r.t. core clock, latency w.r.t. core clock,-->\n']);
        fprintf(fid_xml,['			    <param name="buffer_sizes" value="8, 8, 8, 8"/>	\n']);
        fprintf(fid_xml,['				<!-- all the buffer related are optional -->\n']);
        fprintf(fid_xml,['			    <param name="clockrate" value="1400"/>\n']);
        fprintf(fid_xml,['				<param name="ports" value="1,1,1"/>\n']);
        fprintf(fid_xml,['				<!-- number of r, w, and rw search ports -->\n']);
        fprintf(fid_xml,['				<param name="device_type" value="0"/>\n']);
        fprintf(fid_xml,['				<!-- altough there are multiple access types, \n']);
        fprintf(fid_xml,['				Performance simulator needs to cast them into reads or writes\n']);
        fprintf(fid_xml,['				e.g. the invalidates can be considered as writes -->\n']);
        fprintf(fid_xml,['				<stat name="read_accesses" value="800000"/>\n']);
        fprintf(fid_xml,['				<stat name="write_accesses" value="27276"/>\n']);
        fprintf(fid_xml,['				<stat name="read_misses" value="1632"/>\n']);
        fprintf(fid_xml,['				<stat name="write_misses" value="183"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="20"/>	\n']);
        fprintf(fid_xml,['		</component>\n']);
    end
    for core_number=1:1
        fprintf(fid_xml,['		<component id="system.L2Directory' num2str(core_number-1) '" name="L2Directory' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['				<param name="Directory_type" value="1"/>\n']);
        fprintf(fid_xml,['			    <!--0 cam based shadowed tag. 1 directory cache -->	\n']);
        fprintf(fid_xml,['				<param name="Dir_config" value="1048576,16,16,1,2, 100"/>\n']);
        fprintf(fid_xml,['				<!-- the parameters are capacity,block_width, associativity,bank, throughput w.r.t. core clock, latency w.r.t. core clock,-->\n']);
        fprintf(fid_xml,['			    <param name="buffer_sizes" value="8, 8, 8, 8"/>	\n']);
        fprintf(fid_xml,['				<!-- all the buffer related are optional -->\n']);
        fprintf(fid_xml,['			    <param name="clockrate" value="1400"/>\n']);
        fprintf(fid_xml,['				<param name="ports" value="1,1,1"/>\n']);
        fprintf(fid_xml,['				<!-- number of r, w, and rw search ports -->\n']);
        fprintf(fid_xml,['				<param name="device_type" value="0"/>\n']);
        fprintf(fid_xml,['				<!-- altough there are multiple access types, \n']);
        fprintf(fid_xml,['				Performance simulator needs to cast them into reads or writes\n']);
        fprintf(fid_xml,['				e.g. the invalidates can be considered as writes -->\n']);
        fprintf(fid_xml,['				<stat name="read_accesses" value="58824"/>\n']);
        fprintf(fid_xml,['				<stat name="write_accesses" value="27276"/>\n']);
        fprintf(fid_xml,['				<stat name="read_misses" value="1632"/>\n']);
        fprintf(fid_xml,['				<stat name="write_misses" value="183"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="100"/>	\n']);
        fprintf(fid_xml,['		</component>\n']);
    end
    for core_number=1:num_cores
        fprintf(fid_xml,['		<component id="system.L2cache' num2str(core_number-1) '" name="L2cache' num2str(core_number-1) '">\n']);
        fprintf(fid_xml,['			<!-- all the buffer related are optional -->\n']);
        fprintf(fid_xml,['				<param name="L2_config" value="' num2str(L2_capacity(core_number)) ',' num2str(L2_blockwidth(core_number)) ',' num2str(L2_assoc(core_number)) ',1, ' num2str(L2_latency(core_number)) ',' num2str(L2_latency(core_number)) ', 64,1"/>\n']);
%         fprintf(fid_xml,['				<param name="L2_config" value="' num2str(L2_capacity(core_number)) ',' num2str(64) ',' num2str(L2_assoc(core_number)) ',1, ' num2str(L2_latency(core_number)) ',' num2str(L2_latency(core_number)) ', 64,1"/>\n']);
        fprintf(fid_xml,['				<!-- the parameters are capacity,block_width, associativity, bank, throughput w.r.t. core clock, latency w.r.t. core clock,output_width, cache policy -->\n']);
        fprintf(fid_xml,['				<param name="buffer_sizes" value="' num2str(L2_blockwidth(core_number)) ', ' num2str(L2_blockwidth(core_number)) ', ' num2str(L2_blockwidth(core_number)) ', ' num2str(L2_blockwidth(core_number)) '"/>\n']);
%         fprintf(fid_xml,['				<param name="buffer_sizes" value="' num2str(64) ', ' num2str(64) ', ' num2str(64) ', ' num2str(64) '"/>\n']);
        fprintf(fid_xml,['				<!-- cache controller buffer sizes: miss_buffer_size(MSHR),fill_buffer_size,prefetch_buffer_size,wb_buffer_size-->	\n']);
        fprintf(fid_xml,['				<param name="clockrate" value="' num2str(clk) '"/>\n']);
        fprintf(fid_xml,['				<param name="ports" value="0,0,' num2str(L2_ports(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<!-- number of r, w, and rw ports -->\n']);
        fprintf(fid_xml,['				<param name="device_type" value="0"/>\n']);
        fprintf(fid_xml,['				<stat name="read_accesses" value="' num2str(L2_read_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="write_accesses" value="' num2str(L2_write_accesses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="read_misses" value="' num2str(L2_read_misses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="write_misses" value="' num2str(L2_write_misses(core_number)) '"/>\n']);
        fprintf(fid_xml,['				<stat name="conflicts" value="' num2str(L2_conflicts(core_number)) '"/>	\n']);
        fprintf(fid_xml,['			    <stat name="duty_cycle" value="' num2str(L2_duty(core_number)) '"/>	\n']);
        fprintf(fid_xml,['		</component>\n']);
    end
    fprintf(fid_xml,['\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['<component id="system.L30" name="L30">\n']);
    fprintf(fid_xml,['				<param name="L3_config" value="1048576,64,16,1, 2,100, 64, 1"/>\n']);
    fprintf(fid_xml,['				<!-- the parameters are capacity,block_width, associativity,bank, throughput w.r.t. core clock, latency w.r.t. core clock,-->\n']);
    fprintf(fid_xml,['				<param name="clockrate" value="3500"/>\n']);
    fprintf(fid_xml,['				<param name="ports" value="1,1,1"/>\n']);
    fprintf(fid_xml,['				<!-- number of r, w, and rw ports -->\n']);
    fprintf(fid_xml,['				<param name="device_type" value="0"/>\n']);
    fprintf(fid_xml,['				<param name="buffer_sizes" value="16, 16, 16, 16"/>\n']);
    fprintf(fid_xml,['				<!-- cache controller buffer sizes: miss_buffer_size(MSHR),fill_buffer_size,prefetch_buffer_size,wb_buffer_size-->	\n']);
    fprintf(fid_xml,['				<stat name="read_accesses" value="58824"/>\n']);
    fprintf(fid_xml,['				<stat name="write_accesses" value="27276"/>\n']);
    fprintf(fid_xml,['				<stat name="read_misses" value="1632"/>\n']);
    fprintf(fid_xml,['				<stat name="write_misses" value="183"/>\n']);
    fprintf(fid_xml,['				<stat name="conflicts" value="0"/>	\n']);
    fprintf(fid_xml,['			    <stat name="duty_cycle" value="0.35"/>	\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['		<component id="system.NoC0" name="noc0">\n']);
    fprintf(fid_xml,['			<param name="clockrate" value="' num2str(clk) '"/>\n']);
    fprintf(fid_xml,['			<param name="horizontal_nodes" value="' num2str(num_cores) '"/>\n']);
    fprintf(fid_xml,['			<param name="vertical_nodes" value="1"/>\n']);
    fprintf(fid_xml,['			<param name="has_global_link" value="0"/>\n']);
    fprintf(fid_xml,['			<!-- 1 has global link, 0 does not have global link -->\n']);
    fprintf(fid_xml,['			<param name="link_throughput" value="1"/><!--w.r.t clock -->\n']);
    fprintf(fid_xml,['			<param name="link_latency" value="' num2str(network_latency) '"/><!--w.r.t clock -->\n']);
    fprintf(fid_xml,['			<!-- througput >= latency -->\n']);
    fprintf(fid_xml,['			<!-- Router architecture -->\n']);
    fprintf(fid_xml,['			<param name="input_ports" value="' num2str(router_ports) '"/>\n']);
    fprintf(fid_xml,['			<param name="output_ports" value="' num2str(router_ports) '"/>\n']);
    fprintf(fid_xml,['			<param name="virtual_channel_per_port" value="1"/>\n']);
    fprintf(fid_xml,['			<!-- input buffer; in classic routers only input ports need buffers -->\n']);
    fprintf(fid_xml,['			<param name="flit_bits" value="' num2str(NOC_width) '"/>\n']);
    fprintf(fid_xml,['			<param name="input_buffer_entries_per_vc" value="16"/><!--VCs within the same ports share input buffers whose size is propotional to the number of VCs-->\n']);
    fprintf(fid_xml,['			<param name="chip_coverage" value="1"/>\n']);
    fprintf(fid_xml,['			<!-- When multiple NOC present, one NOC will cover part of the whole chip. chip_coverage <=1 -->\n']);
    fprintf(fid_xml,['			<stat name="total_accesses" value="' num2str(network_accesses) '"/>\n']);
    fprintf(fid_xml,['			<!-- This is the number of total accesses within the whole network not for each router -->\n']);
    fprintf(fid_xml,['		    <stat name="duty_cycle" value="' num2str(network_duty) '"/>\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['		<component id="system.mem" name="mem">\n']);
    fprintf(fid_xml,['			<!-- Main memory property -->\n']);
    fprintf(fid_xml,['			<param name="mem_tech_node" value="' num2str(tech) '"/>\n']);
    fprintf(fid_xml,['			<param name="device_clock" value="' num2str(mem_clk) '"/><!--MHz, this is clock rate of the actual memory device, not the FSB -->\n']);
    fprintf(fid_xml,['			<param name="peak_transfer_rate" value="' num2str(mem_clk*mem_blocksize{1}) '"/><!--MB/S-->\n']);
    fprintf(fid_xml,['			<param name="internal_prefetch_of_DRAM_chip" value="8"/>\n']);
    fprintf(fid_xml,['			<!-- 2 for DDR, 4 for DDR2, 8 for DDR3...-->\n']);
    fprintf(fid_xml,['			<!-- the device clock, peak_transfer_rate, and the internal prefetch decide the DIMM property -->\n']);
    fprintf(fid_xml,['			<!-- above numbers can be easily found from Wikipedia -->\n']);
    fprintf(fid_xml,['			<param name="capacity_per_channel" value="' num2str(mem_size) '"/> <!-- MB -->\n']);
    fprintf(fid_xml,['			<!-- capacity_per_Dram_chip=capacity_per_channel/number_of_dimms/number_ranks/Dram_chips_per_rank\n']);
    fprintf(fid_xml,['			Current McPAT assumes single DIMMs are used.--> 		\n']);
    fprintf(fid_xml,['			<param name="number_ranks" value="2"/>\n']);
    fprintf(fid_xml,['			<param name="num_banks_of_DRAM_chip" value="' num2str(num_banks) '"/>	\n']);
    fprintf(fid_xml,['			<param name="Block_width_of_DRAM_chip" value="' num2str(mem_blocksize{1}) '"/> <!-- B -->\n']);
    fprintf(fid_xml,['			<param name="output_width_of_DRAM_chip" value="' num2str(mem_blocksize{1}) '"/>\n']);
    fprintf(fid_xml,['			<!--number of Dram_chips_per_rank=" 72/output_width_of_DRAM_chip-->\n']);
    fprintf(fid_xml,['			<!--number of Dram_chips_per_rank=" 72/output_width_of_DRAM_chip-->\n']);
    fprintf(fid_xml,['			<param name="page_size_of_DRAM_chip" value="8"/> <!-- 8 or 16 -->\n']);
    fprintf(fid_xml,['			<param name="burstlength_of_DRAM_chip" value="8"/>\n']);
    fprintf(fid_xml,['			<stat name="memory_accesses" value="' num2str(total_mem_accesses) '"/>\n']);
    fprintf(fid_xml,['			<stat name="memory_reads" value="' num2str(total_mem_reads) '"/>\n']);
    fprintf(fid_xml,['			<stat name="memory_writes" value="' num2str(total_mem_writes) '"/>	\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['		<component id="system.mc" name="mc">\n']);
    fprintf(fid_xml,['			<!-- Memeory controllers are for DDR(2,3...) DIMMs -->\n']);
    fprintf(fid_xml,['			<!-- current version of McPAT uses published values for base parameters of memory controller\n']);
    fprintf(fid_xml,['			improvments on MC will be added in later versions. -->\n']);
    fprintf(fid_xml,['			<param name="type" value="0"/> <!-- 1: low power; 0 high performance -->\n']);
    fprintf(fid_xml,['			<param name="mc_clock" value="' num2str(mem_clk) '"/><!--MHz-->\n']);
    fprintf(fid_xml,['			<param name="peak_transfer_rate" value="' num2str(mem_clk*mem_blocksize{1}) '"/><!--MB/S-->\n']);
    fprintf(fid_xml,['			<param name="block_size" value="' num2str(mem_blocksize{1}) '"/><!--(B) the block size of last level cache, which is the unit for one memory burst transfer -->\n']);
    fprintf(fid_xml,['			<param name="number_mcs" value="' num2str(num_MC) '"/>\n']);
    fprintf(fid_xml,['			<!-- current McPAT only supports homogeneous memory controllers -->\n']);
    fprintf(fid_xml,['			<param name="memory_channels_per_mc" value="1"/>\n']);
    fprintf(fid_xml,['			<param name="number_ranks" value="2"/>\n']);
    fprintf(fid_xml,['			<param name="withPHY" value="0"/>\n']);
    fprintf(fid_xml,['			<!-- # of ranks of each channel-->\n']);
    fprintf(fid_xml,['			<param name="req_window_size_per_channel" value="32"/>\n']);
    fprintf(fid_xml,['			<param name="IO_buffer_size_per_channel" value="32"/>\n']);
    fprintf(fid_xml,['			<param name="databus_width" value="128"/>\n']);
    fprintf(fid_xml,['			<param name="addressbus_width" value="51"/>\n']);
    fprintf(fid_xml,['			<!-- McPAT will add the control bus width to the addressbus width automatically -->\n']);
    fprintf(fid_xml,['			<stat name="memory_accesses" value="' num2str(total_mem_accesses) '"/>\n']);
    fprintf(fid_xml,['			<stat name="memory_reads" value="' num2str(total_mem_reads) '"/>\n']);
    fprintf(fid_xml,['			<stat name="memory_writes" value="' num2str(total_mem_writes) '"/>\n']);
    fprintf(fid_xml,['			<!-- McPAT does not track individual mc, instead, it takes the total accesses and calculate\n']);
    fprintf(fid_xml,['			the average power per MC or per channel. This is sufficent for most application. \n']);
    fprintf(fid_xml,['			Further trackdown can be easily added in later versions. --> \n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['		<component id="system.niu" name="niu">\n']);
    fprintf(fid_xml,['			<!-- On chip 10Gb Ethernet NIC, including XAUI Phy and MAC controller  -->\n']);
    fprintf(fid_xml,['			<!-- For a minimum IP packet size of 84B at 10Gb/s, a new packet arrives every 67.2ns. \n']);
    fprintf(fid_xml,['				 the low bound of clock rate of a 10Gb MAC is 150Mhz -->\n']);
    fprintf(fid_xml,['			<param name="type" value="0"/> <!-- 1: low power; 0 high performance -->\n']);
    fprintf(fid_xml,['			<param name="clockrate" value="350"/>\n']);
    fprintf(fid_xml,['			<param name="number_units" value="2"/> <!-- unlike PCIe and memory controllers, each Ethernet controller only have one port -->\n']);
    fprintf(fid_xml,['			<stat name="duty_cycle" value="1.0"/> <!-- achievable max load <= 1.0 -->\n']);
    fprintf(fid_xml,['			<stat name="total_load_perc" value="0.7"/> <!-- ratio of total achived load to total achivable bandwidth  -->\n']);
    fprintf(fid_xml,['			<!-- McPAT does not track individual nic, instead, it takes the total accesses and calculate \n']);
    fprintf(fid_xml,['			the average power per nic or per channel. This is sufficent for most application. -->  			\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['		<component id="system.pcie" name="pcie">\n']);
    fprintf(fid_xml,['			<!-- On chip PCIe controller, including Phy-->\n']);
    fprintf(fid_xml,['			<!-- For a minimum PCIe packet size of 84B at 8Gb/s per lane (PCIe 3.0), a new packet arrives every 84ns. \n']);
    fprintf(fid_xml,['				 the low bound of clock rate of a PCIe per lane logic is 120Mhz -->\n']); 
    fprintf(fid_xml,['			<param name="type" value="0"/> <!-- 1: low power; 0 high performance -->\n']);
    fprintf(fid_xml,['			<param name="withPHY" value="1"/>\n']);
    fprintf(fid_xml,['			<param name="clockrate" value="350"/>\n']);
    fprintf(fid_xml,['			<param name="number_units" value="1"/>\n']);
    fprintf(fid_xml,['			<param name="num_channels" value="8"/> <!-- 2 ,4 ,8 ,16 ,32 -->\n']);
    fprintf(fid_xml,['			<stat name="duty_cycle" value="1.0"/> <!-- achievable max load <= 1.0 -->\n']);
    fprintf(fid_xml,['			<stat name="total_load_perc" value="0.7"/> <!-- Percentage of total achived load to total achivable bandwidth  -->\n']);
    fprintf(fid_xml,['			<!-- McPAT does not track individual pcie controllers, instead, it takes the total accesses and calculate \n']);
    fprintf(fid_xml,['			the average power per pcie controller or per channel. This is sufficent for most application. -->  			\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['		<component id="system.flashc" name="flashc">\n']);
    fprintf(fid_xml,['		    <param name="number_flashcs" value="0"/>\n']);
    fprintf(fid_xml,['			<param name="type" value="1"/> <!-- 1: low power; 0 high performance -->\n']);
    fprintf(fid_xml,['            <param name="withPHY" value="1"/>\n']);
    fprintf(fid_xml,['			<param name="peak_transfer_rate" value="200"/><!--Per controller sustainable reak rate MB/S -->\n']); 
    fprintf(fid_xml,['			<stat name="duty_cycle" value="1.0"/> <!-- achievable max load <= 1.0 -->\n']);
    fprintf(fid_xml,['			<stat name="total_load_perc" value="0.7"/> <!-- Percentage of total achived load to total achivable bandwidth  -->\n']);
    fprintf(fid_xml,['			<!-- McPAT does not track individual flash controller, instead, it takes the total accesses and calculate \n']);
    fprintf(fid_xml,['			the average power per fc or per channel. This is sufficent for most application -->  			\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['<!--**********************************************************************-->\n']);
    fprintf(fid_xml,['\n']);
    fprintf(fid_xml,['		</component>\n']);
    fprintf(fid_xml,['</component>\n']);
    fprintf(fid_xml,['\n']);
    
    fclose(fid_xml);
    
    return

end

