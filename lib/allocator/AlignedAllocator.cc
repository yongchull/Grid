

#include <Grid/GridCore.h>

namespace Grid
{

int PointerCache::victim = 0;
long int PointerCache::TotalAlignedAllocatedBytes = 0;

PointerCache::PointerCacheEntry PointerCache::Entries[PointerCache::Ncache];

PointerCache::PointerInfo PointerCache::Insert(void *ptr, size_t bytes)
{
  PointerInfo pinfo;
  pinfo.address = NULL;
  pinfo.bytes = 0;

  if (bytes < 4096)
  {
    pinfo.address = ptr;
    pinfo.bytes = bytes;
    return pinfo; //delete!
  }

#ifdef GRID_OMP
  assert(omp_in_parallel() == 0);
#endif
  void *ret = NULL;
  int v = -1;

  // Check if all cache entries are valid
  for (uint8_t e =0; e !=Ncache; e++)
    if (Entries[e].valid == false)
    {
      v = e; break;
    }


  // if all valid find a victim
  if (v == -1)
  {
    v = victim;
    victim = (victim + 1) % Ncache;
  }

#if 0
  std::cout << "Victim for cache substitution v="<< v<< std::endl;
  std::cout <<"addr ="<< Entries[v].address << std::endl;
  std::cout <<"bytes="<< Entries[v].bytes << std::endl;
  std::cout <<"valid="<< Entries[v].valid << std::endl;
#endif

  if (Entries[v].valid)
  {
    pinfo.address = Entries[v].address;
    pinfo.bytes = Entries[v].bytes;
    Entries[v].valid   = false;
    Entries[v].address = NULL;
    Entries[v].bytes   = 0;
    //std::cout << "Valid cache entry evicted at address = "<< ret << std::endl;
  }

  Entries[v].address = ptr;
  Entries[v].bytes = bytes;
  Entries[v].valid = true;

#if 0
  std::cout << "New cache entry v="<< v<< std::endl;
  std::cout <<"addr ="<< Entries[v].address << std::endl;
  std::cout <<"bytes="<< Entries[v].bytes << std::endl;
  std::cout <<"valid="<< Entries[v].valid << std::endl;
#endif

  return pinfo;
}

void *PointerCache::Lookup(size_t bytes)
{

  if (bytes < 4096)
    return NULL; // allocate

#ifdef _OPENMP
  assert(omp_in_parallel() == 0);
#endif

  for (uint8_t e = 0; e != Ncache; e++) // faster
  {
    if (Entries[e].valid && (Entries[e].bytes == bytes))
    {
      //std::cout << "Hit a valid entry n." << e << " of bytes " << Entries[e].bytes << " address = " << Entries[e].address << " -- Invalidating" << std::endl;
      Entries[e].valid = 0;
      return Entries[e].address;
    }
  }
  //std::cout << "No valid entries that match" << std::endl;
  return NULL;
}
}
