
inline
static void reverse_bits_of_byte(unsigned char *xx)
{
	unsigned char x = *xx;
	unsigned char t[16] = {
		 0, //  0 = 0000 ~ 0000 =  0
		 8, //  1 = 0001 ~ 1000 =  8
		 4, //  2 = 0010 ~ 0100 =  4
		12, //  3 = 0011 ~ 1100 = 12
		 2, //  4 = 0100 ~ 0010 =  2
		10, //  5 = 0101 ~ 1010 = 10
		 6, //  6 = 0110 ~ 0110 =  6
		14, //  7 = 0111 ~ 1110 = 14
		 1, //  8 = 1000 ~ 0001 =  1
		 9, //  9 = 1001 ~ 1001 =  9
		 5, // 10 = 1010 ~ 0101 =  5
		13, // 11 = 1011 ~ 1101 = 13
		 3, // 12 = 1100 ~ 0011 =  3
		11, // 13 = 1101 ~ 1011 = 11
		 7, // 14 = 1110 ~ 0111 =  7
		15  // 15 = 1111 ~ 1111 = 15
	};

	*xx = (t[x&15] << 4) | t[x>>4];
}

inline
static void switch_2endianness(void *tt, int n)
{
	char *t = tt;
	for (int i = 0; i < n; i++)
	{
		char tmp[2] = {t[0], t[1]};
		t[0] = tmp[1];
		t[1] = tmp[0];
		t += 2;
	}
}

inline
static void switch_3endianness(void *tt, int n)
{
	char *t = tt;
	for (int i = 0; i < n; i++)
	{
		char tmp[3] = {t[0], t[1], t[2]};
		t[0] = tmp[2];
		t[1] = tmp[1];
		t[2] = tmp[0];
		t += 3;
	}
}

inline
static void switch_4endianness(void *tt, int n)
{
	char *t = tt;
	for (int i = 0; i < n; i++)
	{
		char tmp[4] = {t[0], t[1], t[2], t[3]};
		t[0] = tmp[3];
		t[1] = tmp[2];
		t[2] = tmp[1];
		t[3] = tmp[0];
		t += 4;
	}
}

inline
static void switch_8endianness(void *tt, int n)
{
	char *t = tt;
	for (int i = 0; i < n; i++)
	{
		char tmp[8] = {t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7]};
		t[0] = tmp[7];
		t[1] = tmp[6];
		t[2] = tmp[5];
		t[3] = tmp[4];
		t[4] = tmp[3];
		t[5] = tmp[2];
		t[6] = tmp[1];
		t[7] = tmp[0];
		t += 8;
	}
}

