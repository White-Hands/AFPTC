#include<pbc.h>
#include<pbc_test.h>
#define G1 element_init_G1
#define GT element_init_GT
#define Z element_init_Zr
#define t element_t
#define ran element_random
pairing_t pairing;

void init_pairing(){
	char param[2048];
	size_t count = fread(param,1,2048,stdin);
	if (!count) pbc_die("input error");
	pairing_init_set_buf(pairing,param,count);
}


void Zr(element_t n)
{
	Z(n,pairing);
	ran(n);
}

void IG1(element_t n){
	G1(n,pairing);
	ran(n);
}

void IGT(element_t n){
	GT(n,pairing);
	ran(n);
}
