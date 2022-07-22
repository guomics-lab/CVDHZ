# update: 20200303
# author: zhangqiushi


# Delete the row/column who is all NA.
Del_NA_prot = function(matrix,
                       row=T,
                       column=F) {
  temp1 = matrix[!apply(matrix, 1, function(x) sum(!is.na(x)) == 0),]  #apply:for a matrix 1 indicates rows, 2 indicates columns
  temp2 = ifelse(row,
                 list(temp1),
                 list(matrix))
  temp2 = temp2[[1]]

  temp3 = temp2[,!apply(temp2, 2, function(x) sum(!is.na(x)) == 0)]
  temp4 = ifelse(column,
                 list(temp3),
                 list(temp2))
  temp4 = temp4[[1]]
  return(temp4)
}

# show the NA_ratio/minimum/maximum of a matrix.
NA_ratio = function(matrix) {
  sum(is.na(matrix))/ncol(matrix)/nrow(matrix)
}

df_describe = function(matrix) {

  na_num = sum(is.na(matrix))
  na_ratio =  na_num/ncol(matrix)/nrow(matrix)
  na_ratio_in_each_prot = apply(matrix,
                                1, 
                                function(x) {
                                  sum(is.na(x))/ncol(matrix)
                                })
  
  
  v = matrix %>%
    unlist() %>%
    as.vector() %>%
    as.numeric() %>%
    na.omit()
  
  
  cat(paste0("\nProtein: ", nrow(matrix),
             "\nRuns: ", ncol(matrix),
             "\nNA count:", na_num,
             # "\nAll NA proteins count:", 
            #length(names(na_ratio_in_each_prot)[na_ratio_in_each_prot == 1]),
             "\nNA ratio: ", paste0(100*round(na_ratio,
                                              3),
                                   "%"),
             "\nValue range: ",
			 paste(range(v),
				   collapse = "-"),
             "\n"))
  #cat(paste0("All NA proteins:"),names(na_ratio_in_each_prot)[na_ratio_in_each_prot == 1])
  }


# For a matrix filtered by different NA thresholds, see what the number of proteins and NA ratio are.
NA_threshold_table = function(matrix) {
  # NA_ratio = function(matrix) {
  #   sum(is.na(matrix))/nrow(matrix)/ncol(matrix)
  # }
  na_ratio_in_each_prot = apply(matrix, 1, function(x) {
    sum(is.na(x))/ncol(matrix)
  })

  temp = data.frame(sample = names(na_ratio_in_each_prot),
                    na_ratio = na_ratio_in_each_prot,
                    stringsAsFactors = F)

  Table1_na = sapply(10:1, function(x) {
    threshold = x/10

    prot_choose = temp$na_ratio <= threshold
    prot_num = sum(prot_choose)
    na = NA_ratio(matrix[prot_choose,]) %>% round(4)

    return(c(threshold = paste0(10*x,"%"),
             protein_num = prot_num,
             NA_ratio = paste0(100*na,"%")))
  }) %>% t()

  print(Table1_na)
  return(Table1_na)
}

# Filter protein matrix with swissprot fasta.
# The first column of the matrix is entry name.
# fasta <- read_excel("~/uniprot-human-filtered-organism_Human_20367_20200223.xlsx")
# If "dict = t" is selected, two columns of data will be added.(gene and entry)
swissprot_fasta_filter = function(pm,
                                  fasta,
                                  dict = F) {

  fasta <- data.frame(fasta, row.names = fasta$Entry)

  # Delete the protein not in SWISS-prote database.
  pm1 = pm[pm[,1] %in% fasta$Entry,]

  fasta1 = fasta[pm1[,1],]

  fasta1$gene <- sapply(fasta1$Gene.names,
                 function(x) {
                   gene = str_split(x, pattern = " ")[[1]][1]
                   gene = gsub(";","",gene) # swissprot改版，部分gene list用;分隔
				   return(gene)
                   }) %>%
    as.vector()

  paste0("Delete Protein:\n",
         paste(pm[,1][!pm[,1] %in% fasta$Entry],
               collapse = ","),"\n",
         "Deleted protein numbers:\n",
         length(pm[,1][!pm[,1] %in% fasta$Entry])) %>%
    cat()

  fasta1$gene[is.na(fasta1$gene)] = ""

  row.names(pm1) = paste(fasta1$Entry,
                         fasta1$gene,
                         sep = "_")

  pm2 = data.frame(entry = fasta1$Entry,
                   gene = fasta1$gene,
                   pm1[,-1],
                   stringsAsFactors = F)

  pm3 = ifelse(dict,
               list(pm2),
               list(pm2[,-c(1,2)]))[[1]]
  return(pm3)
}


# for txt/tsv/xls etc.
# checkname=F do not work in this function.
read_multiple_txt = function(folder,
                             format) {
  library(stringr)
  library(magrittr)

  p = paste0("*.",
             format,
             "$")
  file_path <- list.files(folder,
                          pattern = p,
                          full.names = T)
  file_name <- list.files(folder,
                          pattern = p,
                          full.names = F)

  temp = list()

  for (i in 1:length(file_path)) {
    if (format %in% c("xlsx","csv")) {

      temp[[i]] = ifelse(format == "xlsx",
                         list(data.frame(read.delim(file_path[i]),
                                         stringsAsFactors = F)),
                         list(data.frame(read_csv(file_path[i]),
                                         stringsAsFactors = F)))[[1]]
    } else {
      temp[[i]] = read.delim(file_path[i]) %>%
        data.frame(stringsAsFactors = F)
    }
  }
  names(temp) = file_name
  return(temp)
}

read_multiple_xlsx = function(folder,
                              format) {

  p = paste0("*.",format, "$")
  file_path <- list.files(folder,
                          pattern = p,
                          full.names = T)
  file_name <- list.files(folder,
                          pattern = p,
                          full.names = F)
  temp = list()
  for (i in 1:length(file_path)) {
    temp[[i]] = read_excel(file_path[i]) %>%
      data.frame(stringsAsFactors = F)
  }
  names(temp) = file_name
  return(temp)
}

# read excel file with multiple sheet
read_multiple_sheet = function(path,
                               colnames = TRUE) {
  library(readxl)
  library(openxlsx)

  sheetName = readxl::excel_sheets(path)
  temp = list()
  for (i in 1:length(sheetName)) {
    temp[[i]] = read.xlsx(path,
                          sheet = i,
                          colNames = colnames)
  }
  names(temp) = sheetName
  return(temp)
}

read_multiple_xlsx_multiple_sheet = function(folder,
                                             format) {

  p = paste0("*.",format, "$")
  file_path <- list.files(folder,
                          pattern = p,
                          full.names = T)
  file_name <- list.files(folder,
                          pattern = p,
                          full.names = F)
  temp = list()
  for (i in 1:length(file_path)) {
    temp[[i]] = read_multiple_sheet(file_path[i])
  }
  names(temp) = file_name
  return(temp)
}

# aggregate by column.
aggregate_column = function(matrix = df,
                            label = label,
                            Fun = "mean",
                            IsMatrixlog2 = T,
                            IsResultLog2 = T) {
  # matrix[matrix == "NaN" | matrix == "Filtered"] = NA

  round(digits = 4) # 为了不让as.numeric命令四舍五入
  df1 = apply(matrix,
              2,
              as.numeric) # the result of spectronaut is character

  df2 = ifelse(IsMatrixlog2,
              list(2^df1),
              list(df1))
  df2 = df2[[1]]

  temp = aggregate(t(df2),
                   by = list(label),
                   FUN = Fun,
                   na.rm = TRUE) %>% t()


  temp1 = temp[-1,] %>%
    matrix(ncol = ncol(temp)) %>%
    apply(2, as.numeric)

  temp2 = ifelse(IsResultLog2,
                 list(log2(temp1)),
                 list(temp1))[[1]]

  temp3 = data.frame(temp2,
                     row.names = row.names(matrix),
                     stringsAsFactors = F)
  names(temp3) = temp[1,]
  return(temp3)
}


# 将matrix按行或列分割为list
as.list_matrix_by_1dim = function(matrix, by_row = F) {

  temp = list()
  for(i in 1:nrow(matrix)) {
    temp[[i]] <- matrix[i,]
  }
  temp1 = list()
  for(i in 1:ncol(matrix)) {
    temp1[[i]] <- matrix[,i]
  }

  result = ifelse(by_row,
                  list(temp),
                  list(temp1))[[1]]
  return(result)
}


# 删除pFind.sepctra里的部分样本，现在集群上布置的有点问题
del_sample_from_pFind.spectra = function(spectra,
                                         del_sampel,
                                         out_file) {
  # del_sample = c("F20190610yingyq_PLA_DDA_11CS3",
  #                "F20190610yingyq_PLA_DDA_11CS4")
  # spectra = "Z:/members/zhangqiushi/PLA/pFind.sepctra/20191219SMED_pFind-Filtered.spectra"
  # out_file = "Z:/members/zhangqiushi/PLA/pFind.sepctra/PLA_del2_202004271611.spectra"
  spec <- read.delim(spectra,
                     stringsAsFactors=FALSE)

  file_name = sapply(spec$File_Name,
                     function(x) {
                       str_split(x,
                                 "\\.")[[1]][1]
                     }) %>% unname() %>%
    unique()
  # length(file_name)

  temp = a413[!file_name %in% del_sample,]

  write.table(temp,
              file = out_file,
              row.names = F,
              col.names = F,
              quote = F,
              sep = "\t")
}


# Volcano
## Note that the label may contain Na values
volcano_plot = function(matrix, # must be data.frame with rowname and colnames. row is the protein or gene.
                        label, # vector or factor.
                        group1,
                        group2,
                        fold_change = 2,
                        fc_need_2power_first = F, # effect the log2fc calculate. T: 2^matrix
                        p_cutoff = 0.05,
                        adjustP = T,
                        output_path = NULL, # Full path
                        adjustP_paired = F,
                        adjustP_var.equal = F) {

  matrix = Del_NA_prot(matrix,T,F)

  l1 = which(label == group1)
  l2 = which(label == group2)
  fc_cutoff = log2(fold_change)

  log2fc_raw = apply(matrix,
                 1,
                 function(x) {
                   log2(mean(na.omit(x[l1])) / mean(na.omit(x[l2])))
                 })
  log2fc_2power = apply(2^matrix,
                        1,
                        function(x) {
                          log2(mean(na.omit(x[l1])) / mean(na.omit(x[l2])))
                        })
  log2fc = ifelse(fc_need_2power_first,
                  list(log2fc_2power),
                  list(log2fc_raw))[[1]]

  p_nodeal <- apply(matrix,
                    1,
                    function(y) {

                      p_try = tryCatch(t.test(y[l1],
                                              y[l2],
                                              paired = adjustP_paired,
                                              var.equal = adjustP_var.equal)$p.value,
                                       error = function(x) NA)
                    })

  p_adjust <- p.adjust(p_nodeal,
                       method="BH")

  p = ifelse(adjustP,
             list(p_adjust),
             list(p_nodeal))[[1]]

  pdf_name1 = paste0("Volcano_",
                     paste0(group1,
                            "_",
                            group2,
                            "_"),
                     nrow(matrix),
                     "prot_",
                     p_cutoff,
                     "pNoAdjust",
                     "_fc",
                     2^fc_cutoff,"_",
                     Sys.Date())
  pdf_name2 = paste0("Volcano_",
                     paste0(group1,
                            "_",
                            group2,
                            "_"),
                     nrow(matrix),
                     "prot_",
                     p_cutoff,
                     "pAdjust",
                     "_fc",
                     2^fc_cutoff,"_",
                     Sys.Date())

  pdf_name = ifelse(adjustP,
                    list(pdf_name2),
                    list(pdf_name1))[[1]]
  ylab = ifelse(adjustP,
                list("-log10 Adjust P value"),
                list("-log10 P value"))[[1]]

  # Let 0 at middle
  n = range(log2fc)
  n = max(abs(n))

  pdf(paste0(output_path,
             pdf_name,
             ".pdf"))
  plot(log2fc,
       -log10(p),
       xlim = c(-n,n),
       col = "#00000033",
       pch = 19,
       main = paste0(group1,
                     "_",
                     group2),
       xlab = "log2 Fold Change",
       ylab = ylab)
  abline(h = -log10(p_cutoff),
         v = c(-fc_cutoff,
               fc_cutoff),
         lty = 2,
         lwd = 1)
  up <- log2fc >= fc_cutoff & p <= p_cutoff
  points(log2fc[up],
         -log10(p[up]), col = 1,
         bg = brewer.pal(9,"YlOrRd")[6],
         pch = 21,
         cex = 2)

  down <- log2fc <= -fc_cutoff & p <= p_cutoff
  points(log2fc[down],
         -log10(p[down]),
         col = 1,
         bg = brewer.pal(11,"RdBu")[9],
         pch = 21,
         cex = 2)
  dev.off()

  # a <<- c(row.names(matrix)[down],
  #         row.names(matrix)[up])
  # eval(parse(text= paste0(i,
  #                         "<<- a")))
  name = data.frame(prot = row.names(matrix),
                    entry = sapply(row.names(matrix),
                                   function(x) {
                                     str_split(x,
                                               "_")[[1]][[1]]
                                   }),
                    gene = sapply(row.names(matrix),
                                  function(x) {
                                    str_split(x,
                                              "_")[[1]][[2]]
                                  }),
                    log2fc = log2fc,
                    p = p_nodeal,
                    p_adjust = p_adjust,
                    # Regulate = NA,
                    stringsAsFactors = F)
  # name$Regulate[up] = "up"
  # name$Regulate[down] = "down"

  name <<- name

  write.xlsx(name,
             paste0(output_path,
                    pdf_name,
                    ".xlsx"),
             row.names = F)
  df_describe(matrix)
  print(table(label))

  cat(paste0("\n",group1,"(",
             sum(label == group1),
             ") / ",group2, "(",
             sum(label == group2),
             "): ",
             "\nUpregulated: ", sum(na.omit(up)),
             "\nDownregulated: ", sum(na.omit(down)),
             "\n"))
}

# 计算n组数据有多少种交集
upset_n = function(n) {
  c_n = sapply(1:n,
               function (x) {
                 choose(n,x)
               })
  sum(c_n)
}

# 将p value转换为*号
p_trans = function(num) {

  symnum.args <- list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                      symbols = c("****", "***", "**", "*", "ns"))

  p.signif <- sapply(num,
                     function(x) {
                       symnum.args$x <- x
                       temp = tryCatch(do.call(stats::symnum,
                                        symnum.args),
                                error = function(x) NA)%>%
                         as.character()
                       return(temp)
                     })

  p.signif[p.signif == "?"] = NA
  return(p.signif)
}


# stop write to txt
# only when the result of sinl.number() is 0. the sink() will stops such diversions.
# The reason is that each sink() in the loop causes sink.number increase
# sapply(1:sink.number(),sink) has something wrong.
# when sink.number() equal to 0. the conection with files still exists.
sink.detach = function() {for (i in 1:sink.number()) sink()}
