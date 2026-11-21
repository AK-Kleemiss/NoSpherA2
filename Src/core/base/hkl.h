#pragma once
#include "nos_types.h"

typedef std::array<int, 3> i3;
typedef std::set<i3> hkl_list;
typedef std::set<i3>::const_iterator hkl_list_it;

typedef std::array<double, 3> d3;
typedef std::array<double, 4> d4;
typedef std::set<d3> hkl_list_d;
typedef std::set<d3>::const_iterator hkl_list_it_d;

struct I3Less {
	bool operator()(const i3& a, const i3& b) const noexcept {
		if (a[0] != b[0]) return a[0] < b[0];
		if (a[1] != b[1]) return a[1] < b[1];
		return a[2] < b[2];
	}
};

//indexed by the major grid point the additional points belong to, this map of RefinePoints has field value and the non-integer index in the grid as second
typedef std::multimap<i3, std::pair<double, d3>, I3Less> Refinepointmap;

template <typename numtype = int>
struct hashFunction
{
	size_t operator()(const std::vector<numtype>& myVector) const
	{
		std::hash<numtype> hasher;
		size_t answer = 0;
		for (numtype i : myVector)
		{
			answer ^= hasher(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
		}
		return answer;
	}
};

template <typename numtype = int>
struct hkl_equal
{
	bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const
	{
		const int size = vec1.size();
		if (size != vec2.size())
			return false;
		int similar = 0;
		for (int i = 0; i < size; i++)
		{
			if (vec1[i] == vec2[i])
				similar++;
			else if (vec1[i] == -vec2[i])
				similar--;
		}
		if (abs(similar) == size)
			return true;
		else
			return false;
	}
};

template <typename numtype = int>
struct hkl_less
{
	bool operator()(const std::vector<numtype>& vec1, const std::vector<numtype>& vec2) const
	{
		if (vec1[0] < vec2[0])
		{
			return true;
		}
		else if (vec1[0] == vec2[0])
		{
			if (vec1[1] < vec2[1])
			{
				return true;
			}
			else if (vec1[1] == vec2[1])
			{
				if (vec1[2] < vec2[2])
				{
					return true;
				}
				else
					return false;
			}
			else
				return false;
		}
		else
			return false;
	}
};
