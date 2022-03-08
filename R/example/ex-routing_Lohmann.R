
file_dir <- system.file("VIC_AnKang/flowdir.tif", package = "hydroTools")
file_fract <- system.file("VIC_AnKang/vic5_fract.asc", package = "hydroTools")

file_soil <- system.file("VIC_AnKang/vic5_soil.csv", package = "hydroTools")
soil <- data.table::fread(file_soil)


info = read_flowdir(file_dir)
info

UH <- Lohmann_UH(file_dir, soil,
    fract = file_fract,
    stn_x = 109.000961,
    stn_y = 32.682355
)
expect_equal(dim(UH$UH), c(96, 17))
