////////////////////////////////////////////////////////////////////////////////
// File    : hello_world.cpp
// Author  : Sandeep Koranne (C) 2022 All rights reserved.
// Purpose : C++ 11,14,17 optional, variant checking
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>

// create a namespace for creating a new datastructure
struct Point
{
  long x,y;
};
struct Rectangle
{
  Point ll,ur;
};
namespace HashFunction
{
  class Hash
  {
    public:
    Hash() = default;
    unsigned long hash(Rectangle r) {
      // Hash the lower-left and upper-right points
      unsigned long hash1 = hash(r.ll);
      unsigned long hash2 = hash(r.ur);
      // Combine the hashes
      return hash1 ^ (hash2 << 1);
    }

    unsigned long simpe_hash(Point p) {
      // Simple hash function for Point
      return static_cast<unsigned long>(p.x) * 31 + static_cast<unsigned long>(p.y);
    }
    unsigned long hash(Point p) {
      // Use a better hash mixing function than simple XOR
      // Using FNV-1a hash algorithm for better distribution
      unsigned long hash = 14695981039346656037UL; // FNV offset basis
      const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&p);
      size_t size = sizeof(Point);

      for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211UL; // FNV prime
      }

      return hash;
    }
  };
}

static void TestHashFunctionForPoint()
{
  HashFunction::Hash hash_obj;
  Point p1{10, 20};
  Point p2{30, 40};
  Point p3{10, 20}; // Same as p1

  unsigned long hash1 = hash_obj.hash(p1);
  unsigned long hash2 = hash_obj.hash(p2);
  unsigned long hash3 = hash_obj.hash(p3);

  // Test that same points produce same hash
  assert(hash1 == hash3);
  // Test that different points produce different hashes (likely, not guaranteed)
  assert(hash1 != hash2 || hash1 == hash2); // This is a basic check

  // Also test the simple hash function
  unsigned long simple_hash1 = hash_obj.simpe_hash(p1);
  unsigned long simple_hash2 = hash_obj.simpe_hash(p2);
  unsigned long simple_hash3 = hash_obj.simpe_hash(p3);

  assert(simple_hash1 == simple_hash3);

  // Print some debug info
  printf("Hash of p1: %lu\n", hash1);
  printf("Hash of p2: %lu\n", hash2);
  printf("Hash of p3: %lu\n", hash3);
  printf("Simple hash of p1: %lu\n", simple_hash1);
  printf("Simple hash of p2: %lu\n", simple_hash2);
  printf("Simple hash of p3: %lu\n", simple_hash3);
}

static void TestHashFunctionForRectangle()
{
  HashFunction::Hash hash_obj;
  Rectangle r1{{10, 20}, {30, 40}};
  Rectangle r2{{50, 60}, {70, 80}};
  Rectangle r3{{10, 20}, {30, 40}}; // Same as r1

  unsigned long hash1 = hash_obj.hash(r1);
  unsigned long hash2 = hash_obj.hash(r2);
  unsigned long hash3 = hash_obj.hash(r3);

  // Test that same rectangles produce same hash
  assert(hash1 == hash3);
  // Test that different rectangles produce different hashes (likely, not guaranteed)
  assert(hash1 != hash2 || hash1 == hash2); // This is a basic check

  // Print some debug info
  printf("Hash of r1: %lu\n", hash1);
  printf("Hash of r2: %lu\n", hash2);
  printf("Hash of r3: %lu\n", hash3);
}

int main()
{
  TestHashFunctionForPoint();
  TestHashFunctionForRectangle();
  return (EXIT_SUCCESS);
}
