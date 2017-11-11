#include <math.h>

/*
 * Construct a normalized moment tensor (i.e. without the scalar M_0)
 *      M = sqrt(2/3)*iso*I + sqrt(1-iso*iso)*[sqrt(1-clvd*clvd)*DC+clvd*CLVD]
 * where
 *      DC_ij = n_i v_j + v_i n_j,
 *      CLVD_ij = (2 N_i N_j - v_i v_j - n_i n_j )/sqrt(3),
 *
 *      1>=iso>=-1,
 *      0.5>=clvd>=-0.5,
 *      from strike, dip, and rake (all in degrees). See A&R P117
 *
 * Source: copy from Prof. Lupei Zhu's gcap/radiats.c
 * Reference: Zhu & Ben-Zion, GJI, 2013; Zhu & Zhou, PCE, 2016
 */
#define DEG2RAD  1.745329252e-2
void nmtensor(float iso, float clvd, float str, float dip, float rake, float tensor[3][3]) {
    float cstr,cdip,crak,sstr,sdip,srak,sstr2,cstr2,sdip2,cdip2,dum,dev;
    float n[3], v[3], N[3];
    dum = 0.8164966*iso;
    tensor[0][0] = tensor[1][1] = tensor[2][2] = dum;
    tensor[0][1] = tensor[0][2] = tensor[1][2] = 0.;
    dev = 1.-iso*iso;
    if (dev>0.) {
        dev = sqrt(dev);
        dum = dev*sqrt(1.-clvd*clvd);
        str  *= DEG2RAD; dip  *= DEG2RAD; rake *= DEG2RAD;
        sstr=sin(str);cstr=cos(str);sstr2=2*sstr*cstr;cstr2=1-2*sstr*sstr;
        sdip=sin(dip);cdip=cos(dip);sdip2=2*sdip*cdip;cdip2=1-2*sdip*sdip;
        crak=cos(rake);srak=sin(rake);
        tensor[0][0] += -dum*(sdip*crak*sstr2+sdip2*srak*sstr*sstr);
        tensor[0][1] +=  dum*(sdip*crak*cstr2+0.5*sdip2*srak*sstr2);
        tensor[0][2] += -dum*(cdip*crak*cstr+cdip2*srak*sstr);
        tensor[1][1] +=  dum*(sdip*crak*sstr2-sdip2*srak*cstr*cstr);
        tensor[1][2] +=  dum*(cdip2*srak*cstr-cdip*crak*sstr);
        tensor[2][2] +=  dum*sdip2*srak;
        if (clvd>0.0001 || clvd<-0.0001) {
            n[0] = -sdip*sstr; n[1] = sdip*cstr; n[2] = -cdip;
            v[0] =  crak*cstr+srak*cdip*sstr; v[1] = crak*sstr-srak*cdip*cstr; v[2] = -srak*sdip;
            N[0] = n[1]*v[2]-n[2]*v[1]; N[1] = n[2]*v[0]-n[0]*v[2]; N[2] = n[0]*v[1]-n[1]*v[0];

            dum = dev*clvd/sqrt(3.);
            tensor[0][0] += dum*(2*N[0]*N[0]-n[0]*n[0]-v[0]*v[0]);
            tensor[0][1] += dum*(2*N[0]*N[1]-n[0]*n[1]-v[0]*v[1]);
            tensor[0][2] += dum*(2*N[0]*N[2]-n[0]*n[2]-v[0]*v[2]);
            tensor[1][1] += dum*(2*N[1]*N[1]-n[1]*n[1]-v[1]*v[1]);
            tensor[1][2] += dum*(2*N[1]*N[2]-n[1]*n[2]-v[1]*v[2]);
            tensor[2][2] += dum*(2*N[2]*N[2]-n[2]*n[2]-v[2]*v[2]);
        }
    }
    tensor[1][0] = tensor[0][1];
    tensor[2][0] = tensor[0][2];
    tensor[2][1] = tensor[1][2];
}
