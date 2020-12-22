#include <iostream>
#include <tommath>

struct coord{
    mp_int x;
    mp_int y;
    mp_int z;
};
void put_in_new_dot_val (coord* c, mp_int* x, mp_int* y, mp_int* z)
{
    mp_init (&(c->x));
    mp_init (&(c->y));
    mp_init (&(c->z));
    mp_copy (x, &(c->x));
    mp_copy (y, &(c->y));
    mp_copy (z, &(c->z));
}
void put_in_dot_val (coord* c, mp_int* x, mp_int* y, mp_int* z)
{
    mp_copy (x, &(c->x));
    mp_copy (y, &(c->y));
    mp_copy (z, &(c->z));
}


void addition(coord* c1, coord* c2, coord* c3, mp_int modulo){
    mp_int t1, t2, t3, t4, t5, t6, t7;
    mp_init_copy (&t1, &(c1->x));
    mp_init_copy (&t2, &(c1->y));
    mp_init_copy (&t3, &(c1->z));
    mp_init_copy (&t4, &(c2->x));
    mp_init_copy (&t5, &(c2->y));
    mp_init_copy (&t6, &(c2->z));
    mp_init (&t7);
    /*T1 = X1
      T2 = Y1
      T3 = Z1
      T4 = X2
      T5 = Y2
      T6 = Z2*/
      //T7 = T1*T6
    mp_mul (&t1, &t6, &t7);
    mp_mod(&t7, &modulo, &t7);
      //T1 = T1*T5
      mp_mul (&t1, &t5, &t1);
      mp_mod(&t1, &modulo, &t1);
      //T5 = T3*T5
    mp_mul (&t3, &t5, &t5);
    mp_mod(&t5, &modulo, &t5);
      //T3 = T3*T4
    mp_mul (&t3, &t4, &t3);
    mp_mod(&t3, &modulo, &t3);
      //T4 = T2*T4
    mp_mul (&t2, &t4, &t4);
    mp_mod(&t4, &modulo, &t4);
      //T2 = T2*T6
    mp_mul (&t2, &t6, &t2);
    mp_mod(&t2, &modulo, &t2);
      //T6 = T2*T7
    mp_mul (&t2, &t7, &t6);
    mp_mod(&t6, &modulo, &t6);
      //T2 = T2*T4
    mp_mul (&t2, &t4, &t2);
    mp_mod(&t2, &modulo, &t2);
      //T4 = T3*T4
    mp_mul (&t3, &t4, &t4);
    mp_mod(&t4, &modulo, &t4);
      //T3 = T3*T5
    mp_mul (&t3, &t5, &t3);
    mp_mod(&t3, &modulo, &t3);
      //T5 = T1*T5
    mp_mul (&t1, &t5, &t5);
    mp_mod(&t5, &modulo, &t5);
      //T1 = T1*T7
    mp_mul (&t1, &t7, &t1);
    mp_mod(&t1, &modulo, &t1);
      //T1 = T1-T4
      mp_sub (&t1, &t4, &t1);
    mp_mod(&t1, &modulo, &t1);

      //T2 = T2-T5
    mp_sub (&t2, &t5, &t2);
    mp_mod(&t2, &modulo, &t2);
      //T3 = T3-T6
    mp_sub (&t3, &t6, &t3);
    mp_mod(&t3, &modulo, &t3);
      //X3 = T2
      put_in_dot_val(c3, &t2, &t1, &t3);
     /* mp_copy (&t2, &(c3->x));
      //Y3 = T1
    mp_copy (&t1, &(c3->y));
      //Z3 = T3
    mp_copy (&t3, &(c3->z));*/

    mp_clear (&t1);
    mp_clear (&t2);
    mp_clear (&t3);
    mp_clear (&t4);
    mp_clear (&t5);
    mp_clear (&t6);
    mp_clear (&t7);


}

coord find (coord* c, mp_int modulo, mp_int numb){
    coord new_c;

    int am_bits = mp_count_bits (&numb);

    struct coord temp_c, zero_c;
    mp_int temp, zero;
    mp_init_multi(&temp, &zero);
    mp_set (&zero, 0);

    mp_int x, y, z;
    mp_init_multi(&x, &y, &z);
    mp_set (&x, -1);
    mp_set (&y, 0);
    mp_set (&z, 1);
    put_in_new_dot_val(&zero_c, &x, &y, &z);
    put_in_new_dot_val(&temp_c, &(c->x), &(c->y), &(c->z));
    for (int i = 0; i<am_bits; ++i){
        mp_set (&temp, 1<<i);
        mp_and(&numb, &temp, &temp);
        if (mp_cmp(&temp, &zero)!=MP_EQ){
            addition(&temp_c, &temp_c, &temp_c, modulo);
            addition(&zero_c, &temp_c, &zero_c, modulo);
        }
        else{
            addition(&temp_c, &zero_c, &temp_c, modulo);
            addition(&zero_c, &zero_c, &zero_c, modulo);
        }
    }
    put_in_new_dot_val(&new_c, &(zero_c.x), &(zero_c.y), &(zero_c.z));
    mp_clear(&x);
    mp_clear(&y);
    mp_clear(&z);
    mp_clear(&temp);
    mp_clear(&zero);
    return new_c;

}


bool check_on_curve (coord* c, mp_int* modulo){
    mp_int x3, y3, z3, d, eq1, eq2, var;
    mp_init_multi (&x3, &y3, &z3, &d, &eq1, &eq2, &var);
    mp_set(&d, 3);
    mp_sqr(&(c->x), &x3);
    mp_mul(&x3, &(c->x), &x3);//xˆ3

    mp_sqr(&(c->y), &y3);
    mp_mul(&y3, &(c->y), &y3);//yˆ3

    mp_sqr(&(c->z), &z3);
    mp_mul(&z3, &(c->z), &z3);//zˆ3

    mp_add(&x3, &y3, &var);
    mp_add(&var, &z3, &eq1);//xˆ3 + yˆ3 + zˆ3

    mp_mul (&(c->x), &(c->y), &var);
    mp_mul (&var, &(c->z), &var);
    mp_mul(&var, &d, &var);
    mp_int koef;
    mp_init (&koef);
    mp_set (&koef, 3);
    mp_mul (&var,  &koef, &eq2); //3*d*x*y*z

    mp_mod(&eq1, modulo, &eq1);
    mp_mod(&eq2, modulo, &eq2);

    bool equal = mp_cmp (&eq1, &eq2) == MP_EQ;
    return equal;


}
bool equal_dots (coord* c1, coord* c2, mp_int modulo){
    mp_int x1, x2, y1, y2, z1, z2;
    mp_init_multi(&x1, &x2, &y1, &y2, &z1, &z2);
    mp_copy (&(c1->x), &x1);
    mp_copy (&(c2->x), &x2);
    mp_copy (&(c1->y), &y1);
    mp_copy (&(c2->y), &y2);
    mp_copy (&(c1->z), &z1);
    mp_copy (&(c2->z), &z2);
    mp_mod (&x1, &modulo, &x1);
    mp_mod (&x2, &modulo, &x2);
    mp_mod (&y1, &modulo, &y1);
    mp_mod (&y2, &modulo, &y2);
    mp_mod (&z1, &modulo, &z1);
    mp_mod (&z2, &modulo, &z2);
    bool eqq = ((mp_cmp(&x1, &x2)==MP_EQ)&&(mp_cmp(&y1, &y2)==MP_EQ)&&(mp_cmp(&z1, &z2)==MP_EQ));
    mp_clear(&x1);
    mp_clear(&x2);
    mp_clear(&y1);
    mp_clear(&y2);
    mp_clear(&z1);
    mp_clear(&z2);
    return eqq;




}

int main() {
    mp_int p, x, y, z, d, m;
    mp_init(&p);
    mp_init(&d);
    mp_init(&x);
    mp_init(&y);
    mp_init(&z);
    mp_init(&m);
    mp_read_radix(&p, "115792089237316195423570985008687907853269984665640564039457584007913111864739",
                  10);
    mp_read_radix(&x, "93528601524582384978654134272795468222405403055717890280271688132874849008326", 10);
    mp_read_radix(&y, "14443324612566128911211262381388707474030458136470034119105598903952521080679", 10);
    mp_set (&d, 3);
    mp_set(&z, 1);
    mp_read_radix(&m, "115792089237316195423570985008687907852907080286716537199505774922796921406320", 10);
    struct coord c;
    put_in_new_dot_val(&c, &x, &y, &z);

    struct coord zero;
    mp_int zero_x, zero_y, zero_z;
    mp_init_multi (&zero_x, &zero_y, &zero_z);
    mp_set(&zero_x, 1);
    mp_set(&zero_y, -1);
    mp_set(&zero_z, 0);
    put_in_new_dot_val(&zero, &zero_x, &zero_y, &zero_z);
    int bits_numb_m = mp_count_bits(&m);
    mp_int rand_num;
    mp_init (&rand_num);
    mp_rand (&rand_num, bits_numb_m-1);

    std::cout<< "Test1\n";
    zero = find (&c, p, rand_num);
    if (check_on_curve(&zero, &p)){
        std::cout<<"Dot does lay on curve\n";
    }
    else { std::cout<<"Dot not does lay on curve\n";}


    struct coord c2;
    mp_set(&zero_x, 1);
    mp_set(&zero_y, -1);
    mp_set(&zero_z, 0);
    put_in_dot_val(&zero, &zero_x, &zero_y, &zero_z);
    put_in_new_dot_val(&c2, &zero_y, &zero_z, &zero_x);
    std::cout<<"Test2\n";
    c2 = find (&c, p, m);

    if (equal_dots(&c2, &zero, p)) {
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}


    std::cout<<"Test3\n";
    mp_copy (&m, &rand_num);
    mp_add (&rand_num, &zero_x, &rand_num);
    c2 = find (&c, p, rand_num);
    if (equal_dots(&c2, &c, rand_num)){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}

    mp_sub (&m, &zero_x, &rand_num);
    c2 = find (&c, p, rand_num);
    mp_set (&zero_x, -1);
    if (equal_dots(&c2, &zero, p)){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}

    std::cout<<"Test4\n";
    mp_int k, k1, k2;
    mp_init_multi(&k1, &k2, &k);
    mp_rand (&k1, 5);
    mp_rand (&k2, 5);
    mp_add (&k1, &k2, &k);
    coord r1 = find (&c, p, k1);
    coord r2 = find (&c, p, k2);
    coord r = find (&c, p, k);
    addition (&r1, &r2, &r1, p);
    if (equal_dots(&r1, &r, p)){
        std::cout<<"OK\n";
    }
    else {std::cout<<"Not OK\n";}
    mp_clear(&k);
    mp_clear(&k1);
    mp_clear(&k2);
    mp_clear(&zero_x);
    mp_clear(&zero_y);
    mp_clear(&zero_z);
    mp_clear(&p);
    mp_clear(&z);
    mp_clear(&x);
    mp_clear(&y);
    mp_clear(&d);
    mp_clear(&m);


    return 0;



    //mp_mulmod
}
