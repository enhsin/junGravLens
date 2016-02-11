
int sign(double x);
double deg2rad(double pha);
void forward_cic(double *cic_in,double *x_in,double *y_in,double bsx,double bsy,int nx,int ny,int np,double *cic_out);
void lanczos_diff_2_tag(double *m1, double *m2, double *m11, double *m12, double *m21, double *m22, double Dcell, int Ncc, int dif_tag);
void xy_rotate(double *x1_in,double *x2_in,int nx1,int nx2,double xc1,double xc2,double pha,double *x1_out,double *x2_out);
void gauss_2d(double *x1,double *x2,int nx1,int nx2,double *par,double *res);
void tophat_2d(double *x1,double *x2, int nx1, int nx2,double *par,double *res);
void lq_nie(double *x1,double *x2,int nx1,int nx2,double *lpar,double *alpha1,double *alpha2);
void find_critical_curve(double *mu,int nx,int ny,double* res);
void tot_lq(double *x1, double *x2,int nx1,int nx2,double *lpar, int npars, double *lpars, int nsubs, double *y1, double *y2);
void refine_critical(double * xi1,double * xi2,int nx1,int nx2,double * lpar,int npars,double * lpars, int nsubs,double * critical,int clen, int nfiner, double * yi1,double *yi2);
void lens_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_lens);
void srcs_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_srcs);
void mmbr_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_edge);
void all_about_lensing(double *xi1,double *xi2,int nx1,int nx2,double * spar, int nspars, double * spars, int nssubs, double * lpar,int nlpars,double * lpars,int nlsubs,double *s_image,double *g_lensimage,double *critical,double *caustic);
