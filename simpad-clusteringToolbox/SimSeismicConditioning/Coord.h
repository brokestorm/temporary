#ifndef _COORD
#define _COORD

struct Coord
{
	int x;
	int y;
	int z;

	Coord()
	{
		x = 0;
		y = 0;
		z = 0;
	}

	Coord(int x, int y, int z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}

	void UpdateValues(int x, int y, int z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	
	

	std::size_t operator() (const Coord &p) const 
	{
		return Pair_CantorPairing(p.x, Pair_CantorPairing(p.y, p.z));
	}

	bool operator==(const Coord &other) const
	{
		return (x == other.x && y == other.y && z == other.z);
	}

	static int Pair_CantorPairing(int x, int y)
	{
		int w = x + y;
		int pairr = w * (w + 1) / 2 + y;
		return pairr;
	}
};
#endif // !_COORD
