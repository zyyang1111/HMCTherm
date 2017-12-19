function [ value ] = reference_multi( struct, struct_name, section, property, id )

    section_idx=find(strcmp(struct_name,section));
    section_idx=section_idx(id);
    property_idx=find(strcmp(struct{section_idx}.names,property));
    value=struct{section_idx}.values{property_idx};

end