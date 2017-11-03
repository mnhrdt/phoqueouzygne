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

