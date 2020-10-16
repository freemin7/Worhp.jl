using CBindingGen

include_dir = try
    ENV["JULIA_WORHP_INCLUDE_PATH"]
catch
    "/usr/include/worhp/"
end

shared_library_so = try
    ENV["JULIA_LIBWORHP_SO"]
catch
    "/usr/lib64/libworhp.so"
end

header_files = joinpath.(include_dir,filter(readdir(include_dir)) do name
    endswith(name,".h") && !startswith(name,"monitor")
end)

converted_header = convert_headers(header_files, args = ["-I", include_dir]) do cursor
    header = CodeLocation(cursor).file

    startswith(header, include_dir) || endswith(header,"stdbool.h") || return false
    println(header)
    return true
end

open("bindings.jl", "w+") do io
	generate(io,shared_library_so => converted_header)
end

txt = read("bindings.jl", String)

open("bindings.jl", "w") do f
	write(f, replace(txt, "InitError" => "InitErr"))
end
