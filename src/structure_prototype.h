struct C_INFDS{
//public:
	int Nx, NxState, Ny, Nbeta, Ntheta, Nb, Nmu, NLambda, Nbpar, Nu, Npar0, totalT, Nsubj, Neta, Nbetax, N;
	arma::mat Z, G,  tspan, H, U1, allDelta, thetatild, sytild, EStild, EItild, Iytild, bAdaptParams, OMEGAb, bacc, Sigmae, dSigmaede, dSigmaede2, Sigmab, dSigmabdb, dSigmabdb2, start, startpars, b, sy, ES, EI, Iy, lowBound, upBound, y0, SigmaEta, mu, dmudparMu, dmudparMu2, Lambda, dLambdparLamb, dLambdaparLambd2, P0, par, trueb, Tfilter, tidx, lens, ICb;
	//arma::Mat<double> Tfilter, tidx, lens, ICb;
	double omega, maxT, delt, gainpara, gainparb, errtrol, errtrol1, gainpara1, gainparb1, setAccept, scaleb;
	int alp, maxIterStage1, MAXITER, KKO, IT, isInfo;
	arma::field<arma::mat> fulldt, timeDiscrete, Deltat, meanY, dXstarAll, dXstarAll2, dXtildthetafAll, dXtildthetafAll2, tobs, Y;
	arma::cube Xtild;
	arma::Mat<double> allT;
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
arma::mat calculateTheta(const int isPar, const arma::mat &y, arma::vec &i, struct C_INFDS &InfDS);
C_INFDS getXtildIC3(const int isPar, const int getDxFlag, const int freeIC, struct C_INFDS &InfDS);
void VDPMeas(arma::vec x, int Ny, int Nx, int NxState, arma::mat InfDS_Lambda, arma::mat mu, arma::mat *yPred, arma::mat *Jy);
//arma::vec ekfContinuous10(int Nsubj, const int N, const int Ny, const int Nx, const int Nb, const int NxState,const arma::mat Lambda, int totalT, const arma::mat Sigmae, const arma::mat Sigmab, const arma::mat mu, const arma::mat b, const arma::imat allT, const arma::cube Xtild, arma::cube Y);

arma::mat dynfunICM(const int isPar, const arma::mat &xin, arma::vec &i, const int t, const int isStart, struct C_INFDS &InfDS);
arma::cube dfdxFreeICM(const int isPar, arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdparFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdx2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdxdpFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdpdxFreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
arma::cube dfdpar2FreeIC(arma::mat &xin, arma::vec &i, int t, int isStart, struct C_INFDS &InfDS);
