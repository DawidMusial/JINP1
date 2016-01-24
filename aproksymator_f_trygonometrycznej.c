#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "makespl.h"
#include "piv_ge_solver.h"




/*
 * Wielomiany bazowe: i - stopien funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */

double
fi( int i, double x)
{	
		
	if ( i == 0 ) {
		return 	1;
	}  else if ( i == 1 ){
		return sin(x);
	}else {
		if(i % 2) return  cos( i / 2 * x) + fi( i-1, x );
		else return sin( ( i - 1 ) / 2 * x) + fi( i-1, x );
	} 
}




/* Pierwsza pochodna wielomianu */
double
dfi( int i, double x)
{
	
	
    if ( i == 0 ) {
		return 	0;
	} else if ( i == 1 ){
		return cos(x);
	}else {
        if(i % 2) return  i / 2 * sin( i / 2 * x) + fi( i-1, x );
		else return ( ( i - 1 ) / 2 ) * cos( ( ( i - 1 ) / 2 ) * x) - fi( i-1, x );
    }

}



/* Druga pochodna wielomianu */
double
d2fi( int i, double x)
{	
	
	if ( i == 0 ) {
               return 	0;
    } else if ( i == 1 ){
				return -sin(x);
	}else {
               if(i % 2) return pow( i / 2, 2.0 ) *   cos( i / 2 * x) - fi( i-1, x );
			   else return pow(( i - 1 ) / 2, 2.0 ) * sin( ( i - 1 ) / 2 * x) - fi( i-1, x );
    }


	

}

/* Ta pochodna wielomianu */
double
d3fi( int i, double x)
{	
        if ( i == 0 ) {
                return 	0;
        } else if ( i == 1 ){
		return -cos(x);
	}else {
               if(i % 2) return pow( (i / 2), 3.0 ) * sin( i / 2 * x) -  fi( i-1, x );
			   else return pow( ( i - 1 ) / 2, 3.0 ) * cos( ( i - 1 ) / 2 * x) + fi( i-1, x );
        }


}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double 		a = x[0];
	double		b = x[pts->n-1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  	char *nbEnv= getenv( "APROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) >= 0 )
		nb = atoi( nbEnv ) + 1 ;

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi( i, x[k]) * fi( j, x[k]) );

			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, nb,  y[k] * fi( j, x[k]));
	}
        write_matrix(eqs, stdout);


	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
        write_matrix(eqs, stdout);
	
        if (alloc_spl(spl, nb) == 0) {
                for (i = 0; i < spl->n; i++) {
                        double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
		
			
                        xx+= 10.0*DBL_EPSILON;  /* zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale */
                        spl->f[i] = 0;
                        spl->f1[i] = 0;
                        spl->f2[i] = 0;
                        spl->f3[i] = 0;
                        for (k = 0; k < nb; k++) {
                                double          ck = get_entry_matrix(eqs, k, nb);
                                spl->f[i]  += ck * fi  (k, xx);
                                spl->f1[i] += ck * dfi (k, xx);
                                spl->f2[i] += ck * d2fi(k, xx);
                                spl->f3[i] += ck * d3fi(k, xx);
                        }
                }
        }



}
