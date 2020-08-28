#include<pbc.h>
#include<pbc_test.h>
#include "20191209.h"
#define mul_element element_mul
#define sub_element element_sub
#define pow_zn element_pow_zn
#define add_element element_add

int main()
{
	int n;
n=100;
	printf("n = %d\n",n);
	//att = keyword
	//
	/*************************************************/
	double start,finish,start_all,finish_all;
	init_pairing();


	//System Setup
	start=pbc_get_time();
	start_all=start;
	t g1;
	IG1(g1);
	element_random(g1);
	t g0;
	IG1(g0);
	element_random(g0);
	t alpha,a;
	Zr(alpha);
	Zr(a);
	t H0,H1,H2;
	Zr(H0);
	IG1(H1);
	Zr(H2);

	t MKABE;
	IG1(MKABE);
	pow_zn(MKABE,g1,alpha);

	t Y1;
	IGT(Y1);
	element_pairing(Y1,g1,g1);
	element_pow_zn(Y1,Y1,alpha);

	t g1_ex_a;
	IG1(g1_ex_a);
	element_pow_zn(g1_ex_a,g1,a);

	t x,h;
	Zr(x);
	IG1(h);
	element_pow_zn(h,g0,x);

	t uka;
	Zr(uka);
	t s;
	Zr(s);
	t rka;
	Zr(rka);
	element_sub(rka,x,uka);
	t ta;
	Zr(ta);

	t K;
	IG1(K);
	element_pow_zn(K,g1_ex_a,ta);
	element_mul(K,K,Y1);

	t L;
	IG1(L);
	element_pow_zn(L,g1,ta);

	t Katt[n+1];
	for (int i=1;i<=n;i++){
		Zr(Katt[i]);
		element_from_hash(Katt[i],Katt[i],10);
	}
	finish= pbc_get_time()-start;
//	printf("Key generation :%.3f ms\n",finish*1000);

	// Task Encryption
	//Content encryption
	start=pbc_get_time();
	t sigma;
	IGT(sigma);

	t lambdai[n+1];
	t Ci[n+1];
	t s_neg;
	Zr(s_neg);
	element_neg(s_neg,s);
	t a_lambdai;
	Zr(a_lambdai);
//	start=pbc_get_time();
	for(int i =1; i<=n;i++){
		Zr(lambdai[i]);
		IG1(Ci[i]);
		element_from_hash(H1,Katt[i],i);
		mul_element(a_lambdai,a,lambdai[i]);
		element_pow2_zn(Ci[i],g1,a_lambdai,H1,s_neg);
	}
//	finish = pbc_get_time();
//	printf("****************\n%f********",finish-start);
	t C;
	IGT(C);
	pow_zn(C,Y1,s);
	mul_element(C,C,sigma);

	t C_pie;
	IG1(C_pie);
	element_from_hash(H2,sigma,10);
	pow_zn(C_pie,g1,H2);

	t C0;
	IG1(C0);
	pow_zn(C0,g1,s);
	finish=pbc_get_time()-start;
	printf("task encryption by customer :%.3fms\n",finish*1000);

	// Keyword encryption
	start=pbc_get_time();
	t ri[n+1];
	t wi1_pie[n+1];
	t wi2_pie[n+1];
	t wi3_pie[n+1];
	t wi1[n+1];
	t wi2[n+1];
	t h_ex_ri;
	IG1(h_ex_ri);
	for(int i=1;i<=n;i++){
		Zr(ri[i]);
		IG1(wi1_pie[i]);
		IG1(wi2_pie[i]);
		pow_zn(wi1_pie[i],g0,ri[i]);
		pow_zn(wi2_pie[i],wi1_pie[i],uka);
		pow_zn(h_ex_ri,h,ri[i]);
		Zr(wi3_pie[i]);
		element_from_hash(wi3_pie[i],h_ex_ri,i);
	}
	finish=pbc_get_time()-start;
	printf("Keyword encryption by customer :%.3fms\n",finish*1000);


	//Keyword re-encryption
	for(int i=1;i<=n;i++){
		IG1(wi1[i]);
		pow_zn(wi1[i],wi1_pie[i],rka);
		mul_element(wi1[i],wi2_pie[i],wi1[i]);
		element_init_same_as(wi2[i],wi3_pie[i]);
	}

	//Interest Encryption
	start=pbc_get_time();
	t tj[n+1];
	t tdj1_pie[n+1];
	t tdj2_pie[n+1];
	t ukb;
	Zr(ukb);
	t ukb_tj;
	Zr(ukb_tj);
	t tj_neq;
	Zr(tj_neq);
	for(int i=1;i<=n;i++){
		Zr(tj[i]);
		IG1(tdj1_pie[i]);
		IG1(tdj2_pie[i]);
		pow_zn(tdj1_pie[i],g0,tj[i]);

		element_neg(tj_neq,tj[i]);
		mul_element(ukb_tj,ukb,tj[i]);
		element_pow2_zn(tdj2_pie[i],h,tj_neq,g0,ukb_tj);
	}
	finish = pbc_get_time()-start;
	printf("Interest Encryption by worker:%.3f ms\n",finish*1000);

	t rkb;
	Zr(rkb);
	t TDj[n+1];
	for (int i=1;i<=n;i++){
		IG1(TDj[i]);
		pow_zn(TDj[i],tdj1_pie[i],rkb);
		mul_element(TDj[i],TDj[i],tdj2_pie[i]);
	}
	//Task Allocation
	start= pbc_get_time();
	t wi1_TDj;
	IG1(wi1_TDj);
	for(int i=1;i<=n;i++){
		element_invert(wi1_TDj,TDj[i]);
		mul_element(wi1_TDj,wi1[i],TDj[i]);
		element_from_hash(H0,wi1_TDj,i);
	}
	finish = pbc_get_time()-start;
	printf("Task Allocation :%.3f ms\n",finish*1000);

	//Ability verification
	t beta;
	Zr(beta);
	t V;
	IG1(V);
	pow_zn(V,g1,beta);

	start = pbc_get_time();
	t shang;
	IGT(shang);
	element_pairing(shang,K,C0);
	mul_element(shang,C,shang);
	t xia1[n+1];
	t xia2[n+1];

	t etaj[n+1];
	t xia;
	IGT(xia);
	element_set1(xia);
	//printf("1\n");
	for(int i=1;i<=n;i++){
		IGT(xia1[i]);
		IGT(xia2[i]);
		element_pairing(xia1[i],L,Ci[i]);
	//printf("2\n");
	// element_printf("C0 %B\n",C0);
	// element_printf("Katt %B\n",K);
		element_pairing(xia2[i],C0,K);
	//printf("3\n");
		mul_element(xia1[i],xia1[i],xia2[i]);
		Zr(etaj[i]);
		pow_zn(xia1[i],xia1[i],etaj[i]);
		mul_element(xia,xia,xia1[i]);
	}
	element_div(shang,shang,xia);
	finish=pbc_get_time()-start;
	printf("Decryption is :%.3f ms\n",finish*1000);

	element_from_hash(H2,sigma,10);
	t P;
	IG1(P);
	pow_zn(P,V,H2);
	pow_zn(P,C_pie,beta);

//	finish_all=pbc_get_time()-start_all;
	//printf("total_time: %.3f\n",finish_all*1000);
	return 0;
}

