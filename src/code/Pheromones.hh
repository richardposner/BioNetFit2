/*============================================================================
// Name        : Pheromones.hh
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/

#ifndef CODE_PHEROMONES_HH_
#define CODE_PHEROMONES_HH_

#include <iostream>
#include <string.h>
#include <stdio.h>

#include <mpi.h> //Raquel classic MPI support

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/interprocess/ipc/message_queue.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/set.hpp>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>

#include <boost/lexical_cast.hpp>



/*
 * Communication class used for communication between Swarm members.
 *
 * We're using boost::IPC for communication between processes when the entire Swarm is run
 * on a single computer. And we're using boost::mpi for communication between processes
 * when the Swarm is running on a cluster.
 *
 * The communicator can be used to share simple messages, such as simulation results and
 * particle states between particles. It can also be used to share complex Swarm objects,
 * such as a Model or Exp object.
 *
 * The class provides generalized communication methods which allow Swarm communication to
 * take place without the Swarm needing to know if we're on a computer or a cluster.
 */

#include "Swarm.hh"

class Swarm;

#define NONBLOCK FALSE
#define BLOCK TRUE

#define ANY_TAG 10
#define GET_RUNNING_PARTICLES 11
#define SEND_RUNNING_PARTICLES 12
#define SIMULATION_END 13
#define SIMULATION_FAIL 14
#define INIT_BREEDING 15
#define DO_BREED 16
#define RECIPROCATE_BREED 16
#define DONE_BREEDING 17
#define NEXT_GENERATION 18
#define GET_PARAMS_FROM_PARTICLE 19
#define FINISHED_WITH_FIT 20
#define SEND_FINAL_PARAMS_TO_PARTICLE 21
#define FIT_FINISHED 22
#define SEND_NUMFLIGHTS_TO_PARTICLE 23
#define MASTER_DIED 24
#define NEW_BOOTSTRAP 25
#define BEGIN_NELDER_MEAD 26
#define END_NELDER_MEAD 27
#define MESSAGE_END 1000

class Pheromones {
public:
	Pheromones();
	~Pheromones();

	void init(Swarm *s);

	struct swarmMessage {
		std::string tag;

		int id;
		int sender;

		std::vector<std::string> message;

		friend class boost::serialization::access;

		template<class Archive>
		void serialize(Archive & ar, const unsigned int version)
		{
			ar & tag;

			ar & id;
			ar & sender;

			ar & message;
		}
	};

	typedef std::unordered_multimap<int, swarmMessage> swarmMsgHolder;
	typedef swarmMsgHolder::iterator swarmMsgHolderIt;

	void sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message, int messageID = -1);
	//int recvMessage(signed int senderID, const int receiverID, int tag, bool block, std::vector<std::vector<std::string>> &messageHolder, bool eraseMessage = true);
	int recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage = true, int messageID = -1);

	int getRank();

	std::vector<std::string> univMessageSender;
	//std::vector<std::vector<std::string>> univMessageReceiver;
	swarmMsgHolder univMessageReceiver;


private:
	friend class Swarm;
	friend class Particle;

	void clearSwarmMessage(swarmMessage& sm);
	std::string serializeSwarmMessage(swarmMessage sm);
	swarmMessage deserializeSwarmMessage(std::string sm);

	Swarm *swarm_;

	// We use pointers for communication objects because we don't necessarily want to initialize
	// the objects in every case

	// MPI Stuff
	boost::mpi::environment * env_;
	boost::mpi::communicator * world_;

	// IPC Stuff

	// Our shared memory object
	//boost::interprocess::managed_shared_memory *segment_;

	// MQ Stuff
	std::vector<boost::interprocess::message_queue *> smq_;

	/*
	typedef boost::interprocess::allocator<char, boost::interprocess::managed_shared_memory::segment_manager>
	CharAllocator;
	typedef boost::interprocess::basic_string<char, std::char_traits<char>, CharAllocator>
	MyShmString;

	typedef boost::interprocess::allocator<MyShmString, boost::interprocess::managed_shared_memory::segment_manager>
	vecAllocator_;
	typedef boost::interprocess::vector<MyShmString,vecAllocator_>
	MyVector_;

	typedef std::pair<const int, MyVector_>
	ValueType_;

	typedef boost::interprocess::allocator<ValueType_, boost::interprocess::managed_shared_memory::segment_manager>
	ShmemAllocator_;
	typedef boost::interprocess::multimap<int, MyVector_, std::less<int>, ShmemAllocator_>
	MyMap_;

	CharAllocator *     charallocator;
	ShmemAllocator_ * alloc_inst;
	vecAllocator_ * vectorallocator;

	MyShmString *vecString_;

	MyMap_ *swarmMap_;
	MyVector_ *swarmVec_;

	*/

	boost::mpi::request recvRequest_;
	boost::mpi::status	recvStatus_;

	//boost::interprocess::interprocess_mutex *mutex_;

	//void putArrayInSHM(std::vector<std::string> theArray);

};

#endif /* CODE_PHEROMONES_HH_ */
