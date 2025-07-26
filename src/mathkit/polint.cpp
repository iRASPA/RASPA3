module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>
#endif

module polint;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

/*
void Interpolation::polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
        int i,m,ns=0;
        double den,dif,dift,ho,hp,w;
        double *c,*d;

        c=(double *)malloc((n)*sizeof(double));
        d=(double *)malloc((n)*sizeof(double));

        dif = std::fabs(x - xa[0]);
        for (i = 0; i < n; i++)
  {
                if ((dift = std::fabs(x - xa[i])) < dif)
    {
                        ns = i;
                        dif = dift;
                }
                c[i] = ya[i];
                d[i] = ya[i];
        }
        *y = ya[ns];
  --ns;
        for (m = 1; m < n; m++)
  {
                for (i = 0; i < n - m; i++)
    {
                        ho = xa[i] - x;
                        hp = xa[i + m] - x;
                        w = c[i + 1] - d[i];
                        if ((den = ho - hp) == 0.0)
      {
                                fprintf(stderr,"Error in routine polint\n");
                                exit(0);
                        }
                        den = w / den;
                        d[i] = hp * den;
                        c[i] = ho * den;
                }
    if (2 * (ns + 1) < (n - m))
    {
      *dy = c[ns + 1];
    }
    else
    {
      *dy = d[ns];
      --ns;
    }
                *y += *dy;

        }
        free(d);
        free(c);
}*/
