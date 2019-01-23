#include<gmp.h>
#include<stdio.h>
#include<stdlib.h>

int main(int argc, char **argv)
{
	printf("-----\n");
	printf("ARBITRARY PRECISION MULTIPLICATIVE INVERSE\n\n");

	// check the arguments
	if(argc != 3)
	{
		printf("usage:\n");
		printf("\t./multiplicative_inverse.out <argument 1> <argument 2>\n");
		printf("argument 1 = the number whose multiplicative inverse is to be found\n");
		printf("argument 2 = the modulus of the multiplicative inverse calculation\n");
		printf("-----\n");
		return 1;
	}

	// read second input argument
	mpz_t m;
	mpz_init(m);
	if(mpz_set_str(m, argv[2], 10))
	{
		printf("%s is not a valid natural number.\n", argv[2]);
		printf("-----\n");
		return 2;
	}

	// check if it is allowed
	if(mpz_cmp_si(m, 1) <= 0)
	{
		mpz_out_str(stdout, 10, m);
		printf(" is not a valid modulus.\n");
		printf("-----\n");
		return 3;
	}

	// read first input argument
	mpz_t s;
	mpz_init(s);
	if(mpz_set_str(s, argv[1], 10))
	{
		printf("%s is not a valid natural number.\n", argv[1]);
		printf("-----\n");
		return 4;
	}

	// show input
	printf("You entered the following.\n");
	printf("number  = ");
	mpz_out_str(stdout, 10, s);
	printf("\n");
	printf("modulus = ");
	mpz_out_str(stdout, 10, m);
	printf("\n\n");

	// if negative, make it positive and congruent to some number in the complete residue set
	printf("Correcting input argument.\n");
	mpz_mod(s, s, m);

	// display corrections, if any
	printf("number  = ");
	mpz_out_str(stdout, 10, s);
	printf("\n");
	printf("modulus = ");
	mpz_out_str(stdout, 10, m);
	printf("\n\n");

	// check if number is zero
	if(mpz_cmp_si(s, 0) == 0)
	{
		printf("0 does not have a multiplicative inverse with respect to any modulus.\n");
		printf("-----\n");
		return 5;
	}

	// start operation assuming that the two arguments are the first two remainders
	mpz_t small; // small remainder is number
	mpz_init(small);
	mpz_set(small, s);
	mpz_t large; // large remainder is modulus
	mpz_init(large);
	mpz_set(large, m);

	// array to store remainders
	mpz_t *r;
	r = malloc(2 * sizeof *r);
	mpz_init(r[0]);
	mpz_set(r[0], m);
	mpz_init(r[1]);
	mpz_set(r[1], s);

	// array to store multipliers
	mpz_t *a;
	a = malloc(0);

	// remainder obtained at each iteration of loop
	mpz_t remainder;
	mpz_init(remainder);
	mpz_set_si(remainder, 1);

	// multiplier obtained at each iteration of loop
	mpz_t multiplier;
	mpz_init(multiplier);

	// count the steps of the loop
	int count;
	count = 0;

	// Euclid's Division Lemma
	printf("Executing Euclid\'s algorithm.\n");
	while(mpz_cmp_ui(remainder, 0))
	{
		// obtain the remainder and multiplier
		mpz_fdiv_qr(multiplier, remainder, large, small);

		// show the output
		printf("#%6d:\t", ++count);
		mpz_out_str(stdout, 10, remainder);
		printf(" = ");
		mpz_out_str(stdout, 10, large);
		printf(" - ");
		mpz_out_str(stdout, 10, multiplier);
		printf(" * ");
		mpz_out_str(stdout, 10, small);
		printf("\n");

		// set up the next iteration
		mpz_set(large, small);
		mpz_set(small, remainder);

		// store remainder in array
		r = realloc(r, (count + 2) * sizeof *r);
		mpz_init(r[count + 1]);
		mpz_set(r[count + 1], remainder);

		// store multiplier in array
		a = realloc(a, count * sizeof *a);
		mpz_init(a[count - 1]);
		mpz_set(a[count - 1], multiplier);
	}
	printf("\n");

	// release the numbers which will not be used again
	mpz_clear(small);
	mpz_clear(large);
	mpz_clear(remainder);
	mpz_clear(multiplier);

	// penultimate remainder is GCD
	printf("greatest common divisor of\n");
	printf("(");
	mpz_out_str(stdout, 10, s);
	printf(", ");
	mpz_out_str(stdout, 10, m);
	printf(")\n");
	printf("=\n");
	mpz_out_str(stdout, 10, r[count]);
	printf("\n\n");

	if(mpz_cmp_ui(r[count], 1))
	{
		printf("The multiplicative inverse of\n");
		mpz_out_str(stdout, 10, s);
		printf("\nwith respect to the modulus\n");
		mpz_out_str(stdout, 10, m);
		printf("\ndoes not exist.\n");
		return 6;
	}

	// now find the inverse
	mpz_t p; // first coefficient
	mpz_init(p);
	mpz_set_ui(p, 1);
	mpz_t q; // second coefficient
	mpz_init(q);
	mpz_neg(q, a[count - 2]);
	mpz_t temp; // to allow swapping and processing
	mpz_init(temp);

	// loop to find inverse
	printf("Calculating multiplicative inverse.\n");
	printf("1");
	while(count-- - 2)
	{
		// get the next coefficients
		mpz_set(temp, p);
		mpz_set(p, q);
		mpz_mul(q, q, a[count - 2]);
		mpz_sub(q, temp, q);

		// show the operations using the coefficients
		printf("\t= ");
		mpz_out_str(stdout, 10, p);
		printf(" * ");
		mpz_out_str(stdout, 10, r[count - 2]);
		printf(" + ");
		mpz_out_str(stdout, 10, q);
		printf(" * ");
		mpz_out_str(stdout, 10, r[count - 1]);
		printf("\n");
	}
	printf("\n");
	mpz_mod(q, q, m); // making answer positive

	// show output
	printf("multiplicative inverse of\n");
	mpz_out_str(stdout, 10, s);
	printf("\nwith respect to the modulus\n");
	mpz_out_str(stdout, 10, m);
	printf("\n=\n");
	mpz_out_str(stdout, 10, q);
	printf("\n");

	// release multiple-precision numbers
	mpz_clear(s);
	mpz_clear(m);
	mpz_clear(*r);
	mpz_clear(*a);
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(temp);

	// release the memory
	free(r);
	free(a);

	printf("-----\n");
	return 0;
}
