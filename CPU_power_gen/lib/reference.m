function [ value ] = reference( struct, struct_name, section, property )

    section_idx=find(strcmp(struct_name,section));
    property_idx=find(strcmp(struct{section_idx}.names,property));
    value=struct{section_idx}.values{property_idx};

end