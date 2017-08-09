nlistlf='/DATA/Ephys/Josh/G_init_filetable_DSP.dat'
ofl='/DATA/Ephys/Josh/G_init_DSPtbl.csv'

dt = read.table(nlistlf,header=TRUE)

df = c()

df$unitID = sort(unique(dt$ID))
Nunits    = length(df$unitID)

### pre-define relevant fields
# get current column names of the table but discard file name
nmlst = names(dt)
nmlst = nmlst[nmlst != 'File' & nmlst != 'unitFile' & nmlst != 'Paradigm' & nmlst != 'NTrials' & nmlst != 'NCorrect'& nmlst != 'ID']

for(nm in nmlst){
    df[[nm]] = rep(NA,Nunits)
}

df$unittype = rep(NA,Nunits)
df$usable   = rep(NA,Nunits)

prdgmlst = unique(dt$Paradigm)
for(prdgm in prdgmlst){
    df[[prdgm]] = rep(NA,Nunits)
}

for(cnt in 1:Nunits){
    
    p = which(dt$ID == df$unitID[cnt])

    for(nm in nmlst){
        df[[nm]][cnt] = as.character(dt[[nm]][p[1]])
    }    

    for(prdgm in prdgmlst){
        pnm = which(dt$ID == df$unitID[cnt] & dt$Paradigm == prdgm)
        if(length(pnm)>1){  df[[prdgm]][cnt] = paste(dt$unitFile[pnm],'|',sep='', collapse="")  }
        if(length(pnm)==1){ df[[prdgm]][cnt] = as.character(dt$unitFile[pnm]) }
     }
}

write.table(df,ofl, row.names =FALSE)

