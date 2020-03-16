function mu = zkViaLanczos(matVecFunc, matVecTransFunc, r, c, k)
  Tk = nonsymmetric_lanczos(matVecFunc, matVecTransFunc, r, c, k);
  
end