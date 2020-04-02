extern int nucas;


float A12(float *x){
  return 0.0;
}

float A11(float *x){
  return 1.0;
}

float A22(float *x){
  return 1.0;
}

float A21(float *x){
  return 0.0;
}

float A00(float *x){
  if(nucas == 3){
    return 1.0;
  }else{
    return 0.0;
  }
}

float BN(float *x){
  return 0.0;
}

float FOMEGA(float *x){
  float val;
  const float PI=3.141592;
  switch(nucas){
   case 1:
     val = 32.*( x[1]*(1-x[1]) + x[0]*(1-x[0]) );
     break;
   case 2:
     val = 2.*PI*PI*sin(PI*x[0])*sin(PI*x[1]);
     break;
   case 3:
     val = (1+2.*PI*PI)*cos(PI*x[0])*cos(PI*x[1]);
     break;
   default:
     printf("*** SOLEX : exemple non prevu. Abandon.\n");
     exit(1);
     break;
  }
  return val;
}

float FN(float *x, int condar){
  const float PI=3.141592;
  float val;

  switch (nucas) {
    case 1 :
      switch (condar) {
        case 1 :
          val = -16.*(1-2*x[1])*(1-x[0])*x[0];
          break;
        case 2 :
          val = 16.*(1-2*x[0])*(1-x[1])*x[1];
          break;
        case 3 :
          val = 16.*(1-2*x[1])*(1-x[0])*x[0];
          break;
        case 4 :
          val = -16.*(1-2*x[0])*(1-x[1])*x[1];
          break;
        default :
          printf("*** FN : numéro de bord non prevu. Abandon.\n");
          exit(1);
          break;
      }
      break;
    case 2 :
      switch (condar) {
        case 1 :
          val = -PI*sin(PI*x[0])*cos(PI*x[1]);
          break;
        case 2 :
          val = PI*sin(PI*x[1])*cos(PI*x[0]);
          break;
        case 3 :
          val = PI*sin(PI*x[0])*cos(PI*x[1]);
          break;
        case 4 :
          val = -PI*sin(PI*x[1])*cos(PI*x[0]);
          break;
        default :
          printf("*** FN : numéro de bord non prevu. Abandon.\n");
          exit(1);
          break;
      }
      break;
    case 3 :
      switch (condar) {
        case 1 :
          val = PI*cos(PI*x[0])*sin(PI*x[1]);
          break;
        case 2 :
          val = -PI*cos(PI*x[1])*sin(PI*x[0]);
          break;
        case 3 :
          val = -PI*cos(PI*x[0])*sin(PI*x[1]);
          break;
        case 4 :
          val = PI*cos(PI*x[1])*sin(PI*x[0]);
          break;
        default :
          printf("*** FN : numéro de bord non prevu. Abandon.\n");
          exit(1);
          break;
      }
      break;
    default :
      printf("*** SOLEX : exemple non prevu. Abandon.\n");
      exit(1);
      break;
  }
  return(val);
}

float UD(float *x){ // correspond à la solution exacte des 3 exemples
  const float PI=3.141592;
  float val;

  switch (nucas) {
    case 1 :
      val=16.*x[0]*x[1]*(1-x[0])*(1-x[1]);
      break;
    case 2 :
      val=sin(PI*x[0])*sin(PI*x[1]);
      break;
    case 3 :
      val=cos(PI*x[0])*cos(PI*x[1]);
      break;
    default :
      printf("*** SOLEX : exemple non prevu. Abandon.\n");
      exit(1);
      break;
  }
  return(val);
}
