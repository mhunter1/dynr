switch(Model,
       {func_address=dynr.funcaddress(file=paste0("./demo/RSODEmodelNoCov.c"),verbose=FALSE,model=themodel)}
       ,{func_address=dynr.funcaddress(file=paste0("./demo/RSODEmodel.c"),verbose=FALSE,model=themodel)}
)

if (exists("res")) {
  rm(res)
}

res <- dynr.run(themodel, data,func_address,tfun)