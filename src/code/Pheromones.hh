/*
 * Pheromones.hh
 *
 *  Created on: Jul 22, 2015
 *      Author: brandon
 */

#ifndef CODE_PHEROMONES_HH_
#define CODE_PHEROMONES_HH_

#include <iostream>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

//#include <boost/interprocess/shared_memory_object.hpp>
//#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/string.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

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

// TODO Streamline tag and ID conversion
#define ANY_TAG -10
#define GET_RUNNING_PARTICLES -11
#define SEND_RUNNING_PARTICLES -12
#define SIMULATION_END -13
#define SIMULATION_FAIL -14
#define MESSAGE_END -1000

class Pheromones {
public:
	Pheromones() {}

	void init(Swarm *s);

	void sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message);
	int recvMessage(signed int senderID, int receiverID, int tag, bool block, std::vector<std::vector<std::string>> &messageHolder, bool eraseMessage = true);

	template<typename T>
	T recvFromAll(T & messageHolder);

	std::vector<std::string> univMessageSender;
	std::vector<std::vector<std::string>> univMessageReceiver;

private:
	Swarm *swarm_;

	// MPI Stuff
	boost::mpi::environment * env_;
	boost::mpi::communicator * world_;

	// IPC Stuff
	boost::interprocess::managed_shared_memory *segment_;

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

	MyShmString *vecString_;

	MyMap_ *swarmMap_;
	MyVector_ *swarmVec_;

	std::vector<int> runningParticles_;

	void putArrayInSHM(std::vector<std::string> theArray);

};

#endif /* CODE_PHEROMONES_HH_ */
