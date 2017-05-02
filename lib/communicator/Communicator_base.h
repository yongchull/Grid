
    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/communicator/Communicator_base.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_COMMUNICATOR_BASE_H
#define GRID_COMMUNICATOR_BASE_H

///////////////////////////////////
// Processor layout information
///////////////////////////////////
#ifdef GRID_COMMS_MPI
#include <mpi.h>
#endif
#ifdef GRID_COMMS_MPI3
#include <mpi.h>
#endif
#ifdef GRID_COMMS_MPI3L
#include <mpi.h>
#endif
#ifdef GRID_COMMS_SHMEM
#include <mpp/shmem.h>
#endif

namespace Grid {

class CartesianCommunicator {
  public:    

  // 65536 ranks per node adequate for now
  // 128MB shared memory for comms enought for 48^4 local vol comms
  // Give external control (command line override?) of this

  static const int      MAXLOG2RANKSPERNODE = 16;            
  static uint64_t MAX_MPI_SHM_BYTES;

  // Communicator should know nothing of the physics grid, only processor grid.
  int              _Nprocessors;     // How many in all
  std::vector<int> _processors;      // Which dimensions get relayed out over processors lanes.
  int              _processor;       // linear processor rank
  std::vector<int> _processor_coor;  // linear processor coordinate
  unsigned long _ndimension;

#if defined (GRID_COMMS_MPI) || defined (GRID_COMMS_MPI3) || defined (GRID_COMMS_MPI3L)
  static MPI_Comm communicator_world;
         MPI_Comm communicator;
  typedef MPI_Request CommsRequest_t;
#else 
  typedef int CommsRequest_t;
#endif

  ////////////////////////////////////////////////////////////////////
  // Helper functionality for SHM Windows common to all other impls
  ////////////////////////////////////////////////////////////////////
  // Longer term; drop this in favour of a master / slave model with 
  // cartesian communicator on a subset of ranks, slave ranks controlled
  // by group leader with data xfer via shared memory
  ////////////////////////////////////////////////////////////////////
#ifdef GRID_COMMS_MPI3

  static int ShmRank;
  static int ShmSize;
  static int GroupRank;
  static int GroupSize;
  static int WorldRank;
  static int WorldSize;

  std::vector<int>  WorldDims;
  std::vector<int>  GroupDims;
  std::vector<int>  ShmDims;
  
  std::vector<int> GroupCoor;
  std::vector<int> ShmCoor;
  std::vector<int> WorldCoor;

  static std::vector<int> GroupRanks; 
  static std::vector<int> MyGroup;
  static int ShmSetup;
  static MPI_Win ShmWindow; 
  static MPI_Comm ShmComm;
  
  std::vector<int>  LexicographicToWorldRank;
  
  static std::vector<void *> ShmCommBufs;

#else 
  static void ShmInitGeneric(void);
  static commVector<uint8_t> ShmBufStorageVector;
#endif 

  /////////////////////////////////
  // Grid information and queries
  // Implemented in Communicator_base.C
  /////////////////////////////////
  static void * ShmCommBuf;

  // Isend/Irecv/Wait, or Sendrecv blocking
  enum CommunicatorPolicy_t { CommunicatorPolicyConcurrent, CommunicatorPolicySequential };
  static CommunicatorPolicy_t CommunicatorPolicy;
  static void SetCommunicatorPolicy(CommunicatorPolicy_t policy ) { CommunicatorPolicy = policy; }

  size_t heap_top;
  size_t heap_bytes;

  void *ShmBufferSelf(void);
  void *ShmBuffer(int rank);
  void *ShmBufferTranslate(int rank,void * local_p);
  void *ShmBufferMalloc(size_t bytes);
  void ShmBufferFreeAll(void) ;
  
  ////////////////////////////////////////////////
  // Must call in Grid startup
  ////////////////////////////////////////////////
  static void Init(int *argc, char ***argv);
  
  ////////////////////////////////////////////////
  // Constructor of any given grid
  ////////////////////////////////////////////////
  CartesianCommunicator(const std::vector<int> &pdimensions_in);
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // Wraps MPI_Cart routines, or implements equivalent on other impls
  ////////////////////////////////////////////////////////////////////////////////////////
  void ShiftedRanks(int dim,int shift,int & source, int & dest);
  int  RankFromProcessorCoor(std::vector<int> &coor);
  void ProcessorCoorFromRank(int rank,std::vector<int> &coor);
  
  int                      IsBoss(void)            ;
  int                      BossRank(void)          ;
  int                      ThisRank(void)          ;
  const std::vector<int> & ThisProcessorCoor(void) ;
  const std::vector<int> & ProcessorGrid(void)     ;
  int                      ProcessorCount(void)    ;
  int                      NodeCount(void)    ;

  ////////////////////////////////////////////////////////////////////////////////
  // very VERY rarely (Log, serial RNG) we need world without a grid
  ////////////////////////////////////////////////////////////////////////////////
  static int  RankWorld(void) ;
  static void BroadcastWorld(int root,void* data, int bytes);
  
  ////////////////////////////////////////////////////////////
  // Reduction
  ////////////////////////////////////////////////////////////
  void GlobalSum(RealF &);
  void GlobalSumVector(RealF *,int N);
  void GlobalSum(RealD &);
  void GlobalSumVector(RealD *,int N);
  void GlobalSum(uint32_t &);
  void GlobalSum(uint64_t &);
  void GlobalSum(ComplexF &c);
  void GlobalSumVector(ComplexF *c,int N);
  void GlobalSum(ComplexD &c);
  void GlobalSumVector(ComplexD *c,int N);
  
  template<class obj> void GlobalSum(obj &o){
    typedef typename obj::scalar_type scalar_type;
    int words = sizeof(obj)/sizeof(scalar_type);
    scalar_type * ptr = (scalar_type *)& o;
    GlobalSumVector(ptr,words);
  }
  
  ////////////////////////////////////////////////////////////
  // Face exchange, buffer swap in translational invariant way
  ////////////////////////////////////////////////////////////
  void SendToRecvFrom(void *xmit,
		      int xmit_to_rank,
		      void *recv,
		      int recv_from_rank,
		      int bytes);
  
  void SendRecvPacket(void *xmit,
		      void *recv,
		      int xmit_to_rank,
		      int recv_from_rank,
		      int bytes);
  
  void SendToRecvFromBegin(std::vector<CommsRequest_t> &list,
			   void *xmit,
			   int xmit_to_rank,
			   void *recv,
			   int recv_from_rank,
			   int bytes);
  
  void SendToRecvFromComplete(std::vector<CommsRequest_t> &waitall);

  double StencilSendToRecvFromBegin(std::vector<CommsRequest_t> &list,
				  void *xmit,
				  int xmit_to_rank,
				  void *recv,
				  int recv_from_rank,
				  int bytes);
  
  void StencilSendToRecvFromComplete(std::vector<CommsRequest_t> &waitall);
  void StencilBarrier(void);

  ////////////////////////////////////////////////////////////
  // Barrier
  ////////////////////////////////////////////////////////////
  void Barrier(void);
  
  ////////////////////////////////////////////////////////////
  // Broadcast a buffer and composite larger
  ////////////////////////////////////////////////////////////
  void Broadcast(int root,void* data, int bytes);
  
  template<class obj> void Broadcast(int root,obj &data)
    {
      Broadcast(root,(void *)&data,sizeof(data));
    };

}; 
}

#endif
