function find_field(criteria::Function)
    if criteria(fields)
        return fields
    else
        for ff in fields.fields
            if find_field_bool(criteria,ff)
                return find_field(criteria,ff)
            else
                error("Cant find field.")
            end
        end
    end
end

function find_field{T<:FieldNode}(criteria::Function,f::T)
    if criteria(f)
        return f
    else
        for ff in fields.fields
            if find_field_bool(criteria,ff)
                return find_field(criteria,ff)
            else
                error("Cant find field.")
            end
        end
    end
end

function find_field{T<:Union{ScalarField,VectorField}}(criteria::Function,f::T)
    if criteria(f)
        return f
    else
        error("Cant find field.")
    end
end

function find_field_bool(criteria::Function)
    find_field_bool(criteria,fields)
end


function find_field_bool{T<:FieldNode}(criteria::Function,f::T)
    return criteria(f)?true:any(map(f->find_field_bool(criteria,f),f.fields))
end

function find_field_bool{T<:Union{ScalarField,VectorField}}(criteria::Function,f::T)
    return criteria(f)
end
