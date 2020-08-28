#include<pbc.h>
#include<pbc_test.h>
#include "20191209.h"
#define mul element_mul

//pairing_t pairing;
//void init_pairing(){
//	char param[2048];
//	size_t count = fread(param,1,2048,stdin);
//	if (!count) pbc_die("input error");
//	pairing_init_set_buf(pairing,param,count);
//}

void file_input(){
	// m 
	FILE *fid;
	fid = fopen("message.txt","rb");

	if(fid==0){
		printf("File open error!\n");
		//	return 0;
	}
	// to get the size of this file
	fseek(fid,0,SEEK_END);		// seek for the end of the file
	int n = (int)ftell(fid);	//return the offset from start to end
	rewind(fid);			// put the point on the start again!
	// read the file
	char pos[n];
	if(pos==NULL){
		printf("error in memoryl");
	}
	fread(pos,sizeof(char),n,fid);
	mpz_t m[n];			// use ASCII code
	for(int i=0;i<n;i++){
		mpz_init_set_ui(m[i],(long)pos[i]);	//convert char to integer
	}
	fclose(fid);
}


int main()
{
int n;
n=100;
		printf("n = %d\n",n);
	//att = keyword
	//
	/*************************************************/
	double start,finish,start_all;
	start=pbc_get_time();
	start_all=start;
	init_pairing();
	/*******************
	 * setup
	 *****************/
	//pairing_t pairing;
	//char param[2048];
	//size_t count = fread(param,1,2048,stdin);
	//if (!count) pbc_die("input error");
	//pairing_init_set_buf(pairing,param,count);



	//System Setup
//	start=pbc_get_time();

	element_t g0,x,y,x_hat;
	element_init_G1(g0,pairing);
	element_init_Zr(x,pairing);
	element_init_Zr(y,pairing);
	element_init_Zr(x_hat,pairing);
	element_random(g0);
	element_random(x);
	element_random(y);
	element_random(x_hat);

	element_t PK0,PK1,PK2;
	G1(PK0,pairing);
	element_init_GT(PK1,pairing);
	element_init_G1(PK2,pairing);
	element_pow_zn(PK0,g0,x);		
	element_pow_zn(PK2,g0,x_hat);		
	element_pairing(PK1,g0,g0);
	element_pow_zn(PK1,PK1,y);
	//t &Y=PK1;	// name PK1 = Y!
	t Y;
	element_init_same_as(Y,PK1);

	element_t ti[2*n];
	element_t A[2*n];
	for(int j=1;j<=n;j++)
	{
		element_init_Zr(ti[j],pairing);
		element_random(ti[j]);
		element_init_G1(A[j],pairing);
		element_pow_zn(A[j],g0,ti[j]);
//		element_printf("Ai1 :%B \n",A[j]);
	}

	finish = pbc_get_time()-start;
//	printf("System Setup: %.3f ms\n",finish*1000);
	// Key Generation
	start= pbc_get_time();
	element_t r,b;
	element_init_Zr(r,pairing);
	element_init_Zr(b,pairing);
	element_random(r);
	element_random(b);

	element_t PKFN,g1;
	//element_init_Zr(yr,pairing);	
	element_init_G1(g1,pairing);
	element_init_GT(PKFN,pairing);
	//element_mul(yr,r,y);
	//element_pairing(PKFN,g0,g0);
	element_pow_zn(PKFN,Y,r);
	element_pow_zn(g1,g0,b);

	t PK_FN;
	IGT(PK_FN);
	t s,s_;
	Zr(s);
	Zr(s_);
	element_invert(s_,s);
	element_pow_zn(PK_FN,PKFN,s_);

		t att[2*n];
	for(int j=1;j<=n;j++){
		Zr(att[j]);
	} 
	t u;
	Zr(u);
	t PKEU;
	IGT(PKEU);
	element_pow_zn(PKEU,Y,u);
	t PK_EU;
	IGT(PK_EU);
	element_pow_zn(PK_EU,PKEU,s_);

	t v,z;
	Zr(v);
	Zr(z);
	t xv;
	Zr(xv);
	element_mul(xv,x,v);


	t Ki3[2*n];
	t H1[2*n];
	t x_add_H1_[2*n];
	t ai[2*n];
	t bi[2*n];
	t Ki1[2*n];
	t Ki2[2*n];
	t bi_div_ti[2*n];
	t a_sub_b_div_t[2*n];	
	for (int j=1; j<=n;j++){
		//t Ki3[j];
		IG1(Ki3[j]);
		element_pow2_zn(Ki3[j],g0,xv,g1,z);
		//t H1[j];
		Zr(H1[j]);
		element_from_hash(H1[j],att[j],10);
		//t x_add_H1_[j];
		Zr(x_add_H1_[j]);
		element_add(x_add_H1_[j],x_hat,H1[j]);
		element_invert(x_add_H1_[j],x_add_H1_[j]);
		element_pow_zn(Ki3[j],Ki3[j],x_add_H1_[j]);

		//t ai[j],bi[j];
		Zr(ai[j]);
		Zr(bi[j]);

		//t Ki1[j],Ki2[j];
		IG1(Ki1[j]);
		//t bi_div_ti[j];
		Zr(bi_div_ti[j]);
		element_div(bi_div_ti[j],bi[j],ti[j]);
		element_pow_zn(Ki1[j],g0,bi_div_ti[j]);

		IG1(Ki2[j]);
		//t a_sub_b_div_t[j];	
		Zr(a_sub_b_div_t[j]);
		element_sub(a_sub_b_div_t[j],ai[j],bi[j]);
		element_div(a_sub_b_div_t[j],a_sub_b_div_t[j],ti[j]);
		element_pow_zn(Ki2[j],g0,a_sub_b_div_t[j]);
	}

	t a;
	element_init_Zr(a,pairing);
	element_set0(a);
	for (int j=1;j<=n;j++){
		element_add(a,a,ai[j]);
	}

	t K0,K1,K2,K3;
	IG1(K0);
	IG1(K1);
	IG1(K2);
	IG1(K3);
	t y_add_xv;
	Zr(y_add_xv);
	element_add(y_add_xv,y,xv);
	element_pow_zn(K0,g0,y_add_xv);

	element_pow2_zn(K1,g0,xv,g1,z);

	element_pow_zn(K2,g0,z);

	t y_sub_a;
	Zr(y_sub_a);
	element_sub(y_sub_a,y,a);
	element_pow_zn(K3,g0,y_sub_a);
	
	finish = pbc_get_time()-start;
//	printf("Key Generation:%.3f ms\n",finish*1000);
	// Task Encryption
	start = pbc_get_time();

	t d;
	element_init_Zr(d,pairing);
	element_set_si(d,n-1);	//d = n-1

//	element_sub(d,n,d);		
	t CT;
	Zr(CT);
	t a_hat[2*n];
	for (int i=0;i<=n-1;i++){
		Zr(a_hat[i]);
	}

	//4.3.1 Content Encryption
	t Ci_hat[2*n];
		t fx;
		element_init_Zr(fx,pairing);
		element_set0(fx);
		t i_exponent;
		element_init_Zr(i_exponent,pairing);
			t x_exponent_i;
			element_init_Zr(x_exponent_i,pairing);
	for (int j = 1; j<=n-3; j++){
		IG1(Ci_hat[j]);
		for( int i = 0; i<=n-2;i++){
			element_set_si(i_exponent,i);
			element_pow_zn(x_exponent_i,x,i_exponent);
			element_mul(x_exponent_i,x_exponent_i,a_hat[i]);
			element_add(fx,fx,x_exponent_i);
			//element_free(x_exponent_i);
		}
		element_pow_zn(Ci_hat[j],g0,H1[j]);
		element_mul(Ci_hat[j],PK2,Ci_hat[j]);
		element_pow_zn(Ci_hat[j],Ci_hat[j],fx);
	//	element_free(fx);
	}

	t h;
	Zr(h);

	t kM;
	Zr(kM);

	t CTM;
	IGT(CTM);
	element_pow_zn(CTM,Y,h);
	element_mul_zn(CTM,CTM,kM);

	t C_pie;
	IG1(C_pie);
	element_pow_zn(C_pie,g0,h);

	t C1_pie;
	IG1(C1_pie);
	t h_add_a_hat;
	Zr(h_add_a_hat);
	element_add(h_add_a_hat,h,a_hat[0]);
	element_pow_zn(C1_pie,g0,h_add_a_hat);

	t C2_pie;
	IG1(C2_pie);
	element_pow_zn(C2_pie,g1,h_add_a_hat);
	
	finish=pbc_get_time()-start;
	printf("task encryption by customer :%.3fms\n",finish*1000);

	// 4.3.2 Keyword Encryption of the Customer
	start = pbc_get_time();
	t ri[2*n];
	t Ii1[2*n];
	t Ii2[2*n];
	t ri_add_H1;
	Zr(ri_add_H1);
	t s_sub_ri;
	Zr(s_sub_ri);
	for (int i=1;i<=n;i++){
		Zr(ri[i]);
		element_add(ri_add_H1,ri[i],H1[i]);
		IG1(Ii1[i]);
		
		element_pow_zn(Ii1[i],A[i],ri_add_H1);
//		element_printf("Ai1 :%B \n",A[i]);

		Zr(Ii2[i]);
		element_sub(s_sub_ri,s,ri[i]);
		element_mul(Ii2[i],s_sub_ri,H1[i]);
//		element_printf("Ii2 :%B \n",Ii2[i]);
	}

	t I0;
	IGT(I0);
	element_pow_zn(I0,Y,s);

	finish=pbc_get_time()-start;
	printf("keyword encryption by customer :%.3fms\n",finish*1000);
	// 4.3.3 Interest Encryption
	t eta;
	Zr(eta);
	
	t T1;
	IG1(T1);
	element_pow_zn(T1,K3,eta);
	t Ti1[2*n];
	for(int j=1;j<=n;j++){
		IG1(Ti1[j]);
		element_pow_zn(Ti1[j],Ki2[j],eta);
//		element_printf("Ti1 :%B \n",Ti1[j]);
	}

	start=pbc_get_time();
	t lambda;
	Zr(lambda);

	t T0,T1_pie;
	Zr(T0);
	IG1(T1_pie);
	element_add(T0,u,lambda);
	t s_mul_lambda;
	Zr(s_mul_lambda);
	element_mul(s_mul_lambda,s,lambda);
	element_pow_zn(T1_pie,T1,s_mul_lambda);

	t Ti1_pie[2*n];
	t Ti2[2*n];
	t lambda_div_H1;
	Zr(lambda_div_H1);
	for (int j=1;j<=n;j++){
		IG1(Ti1_pie[j]);
		element_div(lambda_div_H1,lambda,H1[j]);
		element_pow_zn(Ti1_pie[j],Ti1[j],lambda_div_H1);
		
//		element_printf("Ti1 :%B \n",Ti1[j]);
		IG1(Ti2[j]);
		element_pow_zn(Ti2[j],Ki1[j],lambda_div_H1);
//		element_printf("Ki1 :%B \n",Ki1[j]);
//		element_printf("lambda :%B \n",lambda_div_H1);
//		element_printf("Ti2 :%B \n",Ti2[j]);
	}

	finish = pbc_get_time()-start;
	printf("Interest Encrytion :%.3f ms\n",finish*1000);


	t T0_pie;
	Zr(T0_pie);
	element_add(T0_pie,T0,r);

	t eta_invert;
	Zr(eta_invert);
	element_invert(eta_invert,eta);

	t T1_pie_pie;
	IG1(T1_pie_pie);
	element_pow_zn(T1_pie_pie,T1_pie,eta_invert);

	t Ti1_pie_pie[2*n];
	t Ti2_pie_pie[2*n];
	for (int i = 1;i<=n;i++){
		IG1(Ti2_pie_pie[i]);
		element_set(Ti2_pie_pie[i],Ti2[i]);	//指针
//		element_printf("Ti1_pie_pie_invert :%B \n",Ti2[i]);
	}
	for (int j = 1;j<=n;j++){
		IG1(Ti1_pie_pie[j]);
		element_pow_zn(Ti1_pie_pie[j],Ti1_pie[j],eta_invert);
	}
	// Task Allocation
	t Ii_star[2*n];
	t p1[2*n];
	t p2[2*n];
	t p3[2*n];
	t p4[2*n];
	t lambda_mul_a_sub_b;
	Zr(lambda_mul_a_sub_b);
	t lambda_mul_b;
	Zr(lambda_mul_b);
	t g0_s;
	t g0_lambda_a_b;
	t g0_lambda_b;
	IG1(g0_lambda_a_b);
	IG1(g0_s);
	IG1(g0_lambda_b);
	for (int j=1;j<=n;j++){
		IG1(Ii_star[j]);
		element_pow_zn(Ii_star[j],A[j],Ii2[j]);
		element_mul(Ii_star[j],Ii_star[j],Ii1[j]);

		IGT(p1[j]);
		IGT(p2[j]);
		IGT(p3[j]);
		IGT(p4[j]);
		element_pairing(p1[j],Ii_star[j],Ti1_pie_pie[j]);
		
		element_sub(lambda_mul_a_sub_b,ai[j],bi[j]);
		element_mul(lambda_mul_a_sub_b,lambda_mul_a_sub_b,lambda);
		element_pow_zn(g0_s,g0,s);
		element_pow_zn(g0_lambda_a_b,g0,lambda_mul_a_sub_b);
		element_pairing(p2[j],g0_lambda_a_b,g0_s);

		element_pairing(p3[j],Ii_star[j],Ti2_pie_pie[j]);

		element_mul(lambda_mul_b,lambda,bi[j]);
		element_pow_zn(g0_lambda_b,g0,lambda_mul_b);
		element_pairing(p4[j],g0_s,g0_lambda_b);
	}

	start = pbc_get_time();
	t Ti2_pie_pie_invert;
	IG1(Ti2_pie_pie_invert);
	t Ti1_pie_pie_Ti2_pie_pie;
	IG1(Ti1_pie_pie_Ti2_pie_pie);
	t p5;
	IGT(p5);
	element_set1(p5);
	t p5_temp;
	IGT(p5_temp);
	/*********************/
//	double ceshi;
	for (int j=1; j<=n;j++){
		element_invert(Ti2_pie_pie_invert,Ti2_pie_pie[j]);
//		element_printf("Ti1_pie_pie_invert :%B \n",Ti2_pie_pie[j]);
		element_mul(Ti1_pie_pie_Ti2_pie_pie,Ti1_pie_pie[j],Ti2_pie_pie_invert);
		/****/
//		ceshi = pbc_get_time();
		element_pairing(p5_temp,Ii_star[j],Ti1_pie_pie_Ti2_pie_pie);
//		element_printf("Ti1_pie_pie_Ti2_pie_pie :%B \n",Ti1_pie_pie_Ti2_pie_pie);
//		element_printf("Ii_star :%B \n",Ii_star[j]);
//		element_printf("p5_temp :%B \n",p5_temp);
//		finish = pbc_get_time()-ceshi;
//		printf("ceshi :%f \n",finish);
		element_mul(p5,p5,p5_temp);
	}
	element_pairing(p5_temp,g0,T1_pie_pie);
	element_mul(p5,p5,p5_temp);
	finish = pbc_get_time()-start;
	printf("Task Alloction :%.3f ms\n",finish*1000);

	t p6;
	IGT(p6);
	t I0_T0;
	IGT(I0_T0);
	element_pow_zn(I0_T0,I0,T0_pie);
	element_mul(p6,I0_T0,PK_EU);
	element_mul(p6,p6,PK_FN);
// printf("1\n");
	//Decryption
	t phii[2*n];
	for(int j=1;j<=n-3;j++){
		IGT(phii[j]);
		element_pairing(phii[j],Ki3[j],Ci_hat[j]);

// printf("2\n");
		/********
		 * 这个地方需要改进
		 * ******/
		element_pow_zn(phii[j],phii[j],s);	//肯定不行
	}
// printf("3\n");
	t phi;
	IGT(phi);
	element_set1(phi);
	for(int j=1;j<=n-3;j++){
		element_mul(phi,phii[j],phi);
	}

	start = pbc_get_time();

	t epsilon;
	IGT(epsilon);

	t shangmian1;
	IGT(shangmian1);
	element_pairing(shangmian1,K1,C1_pie);
	t xiamian1;
	IG1(xiamian1);
	element_div(xiamian1,C_pie,C1_pie);
	element_pow_zn(xiamian1,xiamian1,b);
	element_mul(xiamian1,xiamian1,C2_pie);
	t xiamian2;
	IGT(xiamian2);
	element_pairing(xiamian2,K2,xiamian1);
	element_mul(xiamian2,phi,xiamian2);
	element_div(epsilon,shangmian1,xiamian2);
	
	t keyM;
	IGT(keyM);
	t fenzi;
	IGT(fenzi);
	t fenmu;
	IGT(fenmu);
	element_mul(fenzi,CTM,phi);
	element_pairing(fenmu,K0,C_pie);
	element_div(keyM,fenzi,fenmu);

	finish=pbc_get_time()-start;
	printf("Decryption is:%.3f ms\n",finish*1000);
	finish=pbc_get_time()-start_all;
//	printf("total is %.3f ms\n",finish*1000);
	return 0;
}

