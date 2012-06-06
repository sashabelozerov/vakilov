#ifndef __mersenne_h
#define __mersenne_h

/* Mersenne Twister random number generator (Matsumora and Nishimura,
 * 1996).  It is a twisted GFSR(624,397) with a period of 2^19937-1.
 */

#define MT_LEN		624
#define MT_IA		397
#define MT_IB		(MT_LEN - MT_IA)
#define UPPER_MASK	0x80000000
#define LOWER_MASK	0x7FFFFFFF
#define MATRIX_A	0x9908B0DF
#define TWIST(b, i, j)	((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s)	(((s) & 1) * MATRIX_A)

static int mt_index;
static unsigned long mt_buffer[MT_LEN];

static void mt_random_init(void)
{
  int i;
  for (i= 0;  i < MT_LEN;  ++i)
    /* substitute pseudo-random seed for non-reproducible sequence */
    mt_buffer[i]= i * 1073741827;
  mt_index= 0;
}

static unsigned long mt_random(void)
{
  unsigned long *b= mt_buffer;
  int idx= mt_index;
  unsigned long s;
  int i;
  
  if (idx == MT_LEN * sizeof(unsigned long))
    {
      idx= 0;
      i= 0;
      for (;  i < MT_IB;  ++i)
	{
	  s= TWIST(b, i, i + 1);
	  b[i]= b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
	}
      for (; i < MT_LEN-1;  ++i)
	{
	  s= TWIST(b, i, i + 1);
	  b[i]= b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
	}
      s= TWIST(b, MT_LEN - 1, 0);
      b[MT_LEN - 1]= b[MT_IA - 1] ^ (s >> 1) ^ MAGIC(s);
    }
  mt_index= idx + sizeof(unsigned long);
  return *(unsigned long *)((unsigned char *)b + idx);
}

#endif /* __mersenne_h */
