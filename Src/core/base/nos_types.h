#pragma once
#include "../pch.h"

// Pre-definition of classes included later
class WFN;
class cell;
class atom;
class GridManager;
class BasisSet;
class tshb_block;
struct asym_atom;
enum PartitionType { Becke, TFVC, Hirshfeld, RI, MBIS, EMBIS };
enum class MultipoleScheme { TFVC, HIRSHFELD, MBIS, EMBIS, NUCLEAR, MULLIKEN, SANDERSON };
enum class RGBIOrbitalBasis { NAO, ANO };

typedef std::complex<double> cdouble;
typedef std::vector<double> vec;
typedef std::vector<vec> vec2;
typedef std::vector<vec2> vec3;
typedef std::vector<int> ivec;
typedef std::vector<ivec> ivec2;
typedef std::vector<ivec2> ivec3;
typedef std::vector<cdouble> cvec;
typedef std::vector<cvec> cvec2;
typedef std::vector<cvec2> cvec3;
typedef std::vector<std::vector<cvec2>> cvec4;
typedef std::vector<bool> bvec;
typedef std::vector<bvec> bvec2;
typedef std::vector<bvec2> bvec3;
//operator[] fails naming the index instead of reading past the end, so a short line in a truncated file is an error, not UB
struct svec : std::vector<std::string>
{
	using std::vector<std::string>::vector;
	svec(const std::vector<std::string> &v) : std::vector<std::string>(v) {}
	svec(std::vector<std::string> &&v) : std::vector<std::string>(std::move(v)) {}
	const std::string &operator[](size_t i) const
	{
		if (i >= size())
			throw std::out_of_range("field " + std::to_string(i + 1) + " missing, the line has only " + std::to_string(size()));
		return std::vector<std::string>::operator[](i);
	}
	std::string &operator[](size_t i) { return const_cast<std::string &>(std::as_const(*this)[i]); }
};
typedef std::vector<std::filesystem::path> pathvec;
typedef std::chrono::high_resolution_clock::time_point _time_point;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> dMatrix1;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> cMatrix1;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent>> iMatrix1;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent>> bMatrix1;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dMatrix2;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> cMatrix2;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> iMatrix2;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> bMatrix2;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrix3;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrix3;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> iMatrix3;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> bMatrix3;
typedef Kokkos::Experimental::mdarray<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrix4;
typedef Kokkos::Experimental::mdarray<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrix4;
typedef Kokkos::Experimental::mdarray<int, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> iMatrix4;
typedef Kokkos::Experimental::mdarray<bool, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> bMatrix4;

typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent>> dMatrixRef1;
typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> dMatrixRef2;
typedef Kokkos::mdspan<double, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> dMatrixRef3;

typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent>> cMatrixRef1;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent>> cMatrixRef2;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef3;
typedef Kokkos::mdspan<cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> cMatrixRef4;
typedef Kokkos::mdspan<const cdouble, Kokkos::extents<unsigned long long, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent>> ccMatrixRef4;
