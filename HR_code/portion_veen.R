# 安装必要的库
install.packages("VennDiagram")
install.packages("readxl")
install.packages("readr")
install.packages("tools")

# 加载必要的库
library(VennDiagram)
library(readxl)
library(readr)
library(tools)

input_file <-  "F:/Eu_peptido/prepare Horiticulture plant journal/figure4/tissue_veen_MS.csv"
inputfile_ext <- tolower(file_ext(input_file))
if (inputfile_ext == "csv"){
  data <- read.csv(input_file)
}else if (inputfile_ext == "xlsx") {
  data <- read_excel(input_file)
}else {
  stop("不支持的文件类型！")
}

data[] <- lapply(data, as.character)
veen_data <- lapply(data, function(x) na.omit(unique(x))) # 去掉每列中的NA值和重复值，只保留唯一元素
# set.seed(123)
random_colors <- sample(colors(), 4)

venn.plot <- venn.diagram(
  x = venn_data,
  category.names = colnames(data),
  filename = NULL,
  output = TRUE,
  cat.cex = 1.5,
  cex = 1.5,
  # fill = c("red", "green", "blue", "yellow"),
  # cat.col = c("red", "green", "blue", "yellow"),
  fill = random_colors,
  cat.col = random_colors,
  scaled = TRUE,
  print.mode = "raw",
  rotation.degree = 0,
  euler.d = FALSE
)
grid.newpage()
grid.draw(venn.plot)
png("D:/Desktop/venn_plot.png", width = 800, height = 800)
dev.off()