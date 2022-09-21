typedef arma::mat (*dynfunICMPtr)(const int is_meanb, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS);
typedef arma::cube (*dfdxFreeICMPtr)(const int is_meanb, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
typedef arma::cube (*cubePtr)(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
typedef void (*setParsFreeICwbPtr)(C_INFDS &InfDS);
typedef arma::mat (*_rcpparma_hello_worldPtr)(void);
typedef void (*setParsFreeICwbPtr)(C_INFDS &InfDS);
typedef arma::mat (*getCovarianceMatrixPtr)(C_INFDS &InfDS);

struct C_FUNC_POINTER {
	dynfunICMPtr dynfunICM;
	dfdxFreeICMPtr dfdxFreeICM;
	cubePtr dfdparFreeIC, dfdx2FreeIC, dfdxdpFreeIC, dfdpdxFreeIC, dfdpar2FreeIC;
	setParsFreeICwbPtr setParsFreeICwb;
	getCovarianceMatrixPtr getCovarianceMatrix;
	_rcpparma_hello_worldPtr test;
};

struct C_INFDS{
//public:
	int Nx, NxState, Ny, Nbeta, Ntheta, Nb, Nmu, NLambda, Nsigmae, Nbpar, Nu, Npar0, totalT, Nsubj, Neta, Nbetax, N, MAXGIB;
	arma::mat Z, G,  tspan, H, U1, allDelta, thetatild, sytild, EStild, EItild, Iytild, bAdaptParams, OMEGAb, bacc, Sigmae, dSigmaede, dSigmaede2, Sigmab, dSigmabdb, dSigmabdb2, start, startpars, b, meanb, useb, sy, ES, EI, Iy, lowBound, upBound, y0, SigmaEta, mu, dmudparMu, dmudparMu2, Lambda, dLambdparLamb, dLambdparLamb2, P0, par, trueb, Tfilter, tidx, lens, ICb;
	//arma::Mat<double> Tfilter, tidx, lens, ICb;
	double omega, maxT, delt, gainpara, gainparb, errtrol, errtrol1, gainpara1, gainparb1, setAccept, scaleb;
	int alp, maxIterStage1, MAXITER, KKO, IT, isInfo, setScaleb;
	arma::field<arma::mat> fulldt, timeDiscrete, Deltat, meanY, dXstarAll, dXstarAll2, dXtildthetafAll, dXtildthetafAll2, dXtildthetafAll_meanb, dXtildthetafAll2_meanb, tobs, Y; 
	arma::field<arma::cube> dXtildAll, d2XtildAll; // note need to modify dynrArmadillo
	arma::cube Xtild, dXtild, d2Xtild, Xtild_meanb;
	arma::Mat<double> allT;
	struct C_FUNC_POINTER fp;
};

struct C_INFDS0{
//public:
	double delt;
	arma::mat trueb;
	arma::field<arma::mat> trueX;
};


/*function prototypes*/
double round2(double x, double y);
double quantile(arma::vec X, double p);
arma::vec randsample(arma::vec &source, int length, int size);
arma::vec span_vec(int start, int end, int step);
arma::mat rowProjection(arma::mat data, arma::mat index);
arma::mat colProjection(arma::mat data, arma::mat index);
arma::mat modProjection(arma::mat data, int divisor, int remainder);
int maxElement(arma::imat integer_matrix);
arma::mat calculateTheta(const int is_meanb, const arma::mat &y, arma::vec &i, struct C_INFDS &InfDS);
C_INFDS getXtildIC3(const int is_meanb, const int getDxFlag, const int freeIC, struct C_INFDS &InfDS);
void meas(arma::vec x, int Ny, int Nx, int NxState, arma::mat InfDS_Lambda, arma::mat mu, arma::mat *yPred, arma::mat *Jy);
//arma::vec ekfContinuous10(int Nsubj, const int N, const int Ny, const int Nx, const int Nb, const int NxState,const arma::mat Lambda, int totalT, const arma::mat Sigmae, const arma::mat Sigmab, const arma::mat mu, const arma::mat b, const arma::imat allT, const arma::cube Xtild, arma::cube Y);

arma::mat dynfunICM(const int is_meanb, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS);
arma::cube dfdxFreeICM(const int is_meanb, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdx2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);

/*
//data structure that contains paramters that needs to be returned to dynr/R.
struct C_OUTPUT{
//public:
  int convFlag, nIterStage1, nIterStage2;
  double ss, avebAccept;
  double *Iytild, *thetatild; //column-major to keep the matrix InfDS.Iytild
};
*/

struct C_OUTPUT{
//public:
  int convFlag, nIterStage1, nIterStage2;
  double ss, avebAccept;
  mat Iytild; 
  vec thetatild; //column-major to keep the matrix InfDS.Iytild
  mat b;
  arma::cube Xtild;
  arma::field<arma::mat> dXtildthetafAll;
  
};


