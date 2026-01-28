#cd(raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator")

using .Masslist_Pro
Masslist_Pro.main()



# # #Masslist debug code
## ------------------------------------------------ debug 1: delete all depulications ------------------------------------------------
using CSV
using DataFrames

input_filename = raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\Test_Tof_data\results\masslist_text_v8_debug.csv"
output_filename = raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\Test_Tof_data\results\masslist_text_v8_debug2.csv"

header_line_count = 7
all_lines = readlines(input_filename)
header_block = all_lines[1:header_line_count]
data_as_string = join(all_lines[header_line_count+1:end], "\n")
df = CSV.read(IOBuffer(data_as_string), DataFrame;
              delim='\t',
              header=names(CSV.read(IOBuffer(header_block[end]), DataFrame, delim='\t')),
              ignoreemptylines=true)
row_hashes = [hash(Tuple(row)) for row in eachrow(df)]
hash_counts = Dict()
for h in row_hashes
    hash_counts[h] = get(hash_counts, h, 0) + 1
end
unique_row_hashes = Set(k for (k, v) in hash_counts if v == 1)
filtered_df = df[in.(row_hashes, Ref(unique_row_hashes)), :]
open(output_filename, "w") do file
    for line in header_block
        write(file, line * "\n")
    end
    for row in eachrow(filtered_df)
        row_as_strings = [ismissing(item) ? "" : string(item) for item in row]
        line_to_write = join(row_as_strings, "\t")
        write(file, line_to_write * "\n")
    end
end

println("\nProcessing finished successfully!")



## ------------------------------------------------ debug 1: delete depulication only ------------------------------------------------
# using CSV
# using DataFrames
# input_filename = raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\Test_Tof_data\results\masslist_text_v8.csv"
# output_filename = raw"D:\BaiduSyncdisk\zhensen_working\Work in Europe\Tools\Masslist_Creator\Test_Tof_data\results\masslist_text_v8_debug.csv"

# header_line_count = 7
# all_lines = readlines(input_filename)
# header_block = all_lines[1:header_line_count]
# data_as_string = join(all_lines[header_line_count+1:end], "\n")
# df = CSV.read(IOBuffer(data_as_string), DataFrame;
#               delim='\t',
#               header=names(CSV.read(IOBuffer(header_block[end]), DataFrame, delim='\t')),
#               ignoreemptylines=true)
# deduplicated_df = unique(df, :Mass)
# open(output_filename, "w") do file
#     for line in header_block
#         write(file, line * "\n")
#     end
#     for row in eachrow(deduplicated_df)
#         row_as_strings = [ismissing(item) ? "" : string(item) for item in row]
#         line_to_write = join(row_as_strings, "\t")
#         write(file, line_to_write * "\n")
#     end
# end
# println("\nProcessing finished successfully!")

