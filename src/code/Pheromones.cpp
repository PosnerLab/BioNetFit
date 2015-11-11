/*
 * Pheromones.cpp
 *
 *  Created on: Jul 23, 2015
 *      Author: brandon
 */

#include "Pheromones.hh"

namespace mpi = boost::mpi;
using namespace boost::interprocess;
using boost::lexical_cast;


// TODO: Need lots of error checking here...

void Pheromones::init(Swarm *s) {
	swarm_ = s;

	// Using MPI
	if (swarm_->options_.useCluster) {
		// Set up our MPI environment and communicator
		mpi::environment *env_ = new mpi::environment();
		mpi::communicator *world_ = new mpi::communicator();
	}
	// Using IPC
	else {
		if (s->getIsMaster()) {
			struct shm_remove
			{
				shm_remove() {
					shared_memory_object::remove("Swarm");
				}
				~shm_remove(){
					//shared_memory_object::remove("Swarm");
					//TODO: Destructor gets called when we go out of scope. Let's call it somewhere else..
				}
			} remover;

			managed_shared_memory *segment = new managed_shared_memory(create_only, "Swarm", 10000000);
			segment_ = segment;
		}
		else {
			managed_shared_memory *segment = new managed_shared_memory(open_only, "Swarm");
			segment_ = segment;
		}

		CharAllocator     charallocator  (segment_->get_segment_manager());
		ShmemAllocator_ alloc_inst (segment_->get_segment_manager());
		vecAllocator_ vectorallocator (segment_->get_segment_manager());


		if (s->getIsMaster()) {
			MyShmString *mystring = new MyShmString(charallocator);
			MyMap_ *mymap = segment_->construct<MyMap_>("SwarmMap")(std::less<int>(),alloc_inst);
			MyVector_ *myshmvec = segment_->construct<MyVector_>("myvec")(vectorallocator);
			vecString_ = mystring;
			swarmMap_ = mymap;
			swarmVec_ = myshmvec;
		}
		else {
			MyShmString *mystring = new MyShmString(charallocator);
			MyMap_ *mymap = segment_->find<MyMap_>("SwarmMap").first;
			MyVector_ *myshmvec = segment_->find<MyVector_>("myvec").first;
			vecString_ = mystring;
			swarmMap_ = mymap;
			swarmVec_ = myshmvec;
		}
	}
}

void Pheromones::sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message) {
	// Using MPI
	if (swarm_->options_.useCluster) {
		std::vector<int> receivers;

		// Sending to the entire swarm
		if (receiverID == -1) {
			// First we need to get a list of all running particles from the master
			world_->send(0, GET_RUNNING_PARTICLES, "");
			world_->recv(0, SEND_RUNNING_PARTICLES, runningParticles_);
			receivers = runningParticles_;
		}
		else {
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		// Loop through receivers and perform the send operation
		for (std::vector<int>::iterator i = receivers.begin(); i != receivers.end(); ++i) {
			// Blocking send
			if (block) {

			}
			// Non-blocking send
			else {
				world_->isend(receiverID, tag, message);
			}
		}
	}
	// Using IPC
	else {
		std::vector<int> receivers;

		// Sending to entire swarm
		if (receiverID == -1) {
			// First we need a list of running particles, so let's construct a request
			// and get the list from the master

			*vecString_ = std::to_string(GET_RUNNING_PARTICLES).c_str(); // Set tag
			swarmVec_->push_back(*vecString_); // Push tag
			*vecString_ = std::to_string(rand()).c_str(); // Set unique ID
			swarmVec_->push_back(*vecString_); // Store unique ID
			*vecString_ = std::to_string(senderID).c_str(); // Set sender
			swarmVec_->push_back(*vecString_); // Push sender
			*vecString_ = std::to_string(0).c_str(); // Set message size of 0
			swarmVec_->push_back(*vecString_); // Store message size

			swarmMap_->insert(std::pair<const int,MyVector_>(0,*swarmVec_)); // Insert vector to first position (master slot) of map
			swarmVec_->clear(); // Clear our message vector

			// Receive a message from master containing list of running particles
			std::vector<std::vector<std::string>> messageHolder;
			recvMessage(0,senderID,SEND_RUNNING_PARTICLES,true,messageHolder);
			//std::cout << "got list of running particles. sending message" << std::endl;

			// [RECEIVER]
			// 		[TAG]
			//		[ID]
			// 		[SENDER]
			//		[MESSAGE_SIZE]
			//		[ARR 0]
			//		[ARR 1]
			//		...

			// Add all running particles to list of message receivers
			for (std::vector<std::string>::iterator r = messageHolder[0].begin(); r != messageHolder[0].end(); ++r) {
				receivers.push_back( atoi(r->c_str()) );
			}
		}
		else {
			// If we're not sending to swarm, add receiverID to list of message receivers
			// In this case, the list will only contain 1 item
			receivers.push_back(receiverID);
		}

		unsigned long int messageID = rand();

		// Add tag, sender, message size, and message(s) if applicable
		for(std::vector<int>::iterator r = receivers.begin(); r != receivers.end(); ++r) {
			int numMsgElements = message.size();

			*vecString_ = std::to_string(tag).c_str(); // Set tag
			swarmVec_->push_back(*vecString_); // Store tag

			*vecString_ = std::to_string(messageID).c_str(); // Set unique ID
			swarmVec_->push_back(*vecString_); // Store unique ID

			*vecString_ = std::to_string(senderID).c_str(); // Set sender
			swarmVec_->push_back(*vecString_); // Store sender

			*vecString_ = std::to_string(numMsgElements).c_str(); // Set message size
			swarmVec_->push_back(*vecString_); // Store message size

			// Insert messages in vector
			for(std::vector<std::string>::iterator m = message.begin(); m != message.end(); ++m) {
				*vecString_ = m->c_str();
				swarmVec_->push_back(*vecString_);
			}

			//std::cout << "storing vector in swarm map to receiver: " << *r << " with tag: " << tag << std::endl;
			// Insert vector to map
			swarmMap_->insert(std::pair<const int,MyVector_>(*r,*swarmVec_));
			string tests = lexical_cast<string>(*(swarmMap_->find(*r)->second.begin()));
			//std::cout << "swarmmap is now: " << swarmMap_->size() << " " <<  tests << std::endl;
			swarmVec_->clear();
		}
		//std::cout << "done sending " << std::endl;

		// Add code to check whether message was seen or not. Only move on once it's gone
		if (block) {
			bool foundMessage = true;
			while (foundMessage) {
				usleep(250000);

				if (recvMessage(senderID, receiverID, tag, false, univMessageReceiver, false)) {
					foundMessage = true;
				}
				else {
					foundMessage = false;
				}
			}
		}
	}
	//std::cout << "really done sending " << std::endl;
}

int Pheromones::recvMessage(signed int senderID, int receiverID, int tag, bool block, std::vector<std::vector<std::string>> &messageHolder, bool eraseMessage) {
	int numMessages = 0;
	if (swarm_->options_.useCluster) {

	}
	else {
		while (1) {
			//std::cout << "sm size: " << swarmMap_->size() << std::endl;
			// TODO: Do we really need to iterate here? Can't just look up by map key??
			// Maybe not.  Think we get a segfault when using map find(). But why?
			for (MyMap_::iterator s = swarmMap_->begin(); s != swarmMap_->end();) {
				//std::cout << receiverID << "looking at " << s->first << std::endl;
				if (s->first == receiverID) {
					//std::cout << "first: " << s->first << " second: " << *(s->second.begin()) << std::endl;
					int currTag;
					int currSender;
					for (MyVector_::iterator v = s->second.begin(); v != s->second.end(); ) {
						//std::cout << receiverID << ": " << *v << std::endl;
						std::string theString = lexical_cast<std::string>(*v);
						// TODO: We have to lexical cast to string, THEN convert to int because a lexical cast to int will throw an exception if the input can't be easily converted to an integer.  Is there a more elegant way to do this?
						if (!theString.empty() && stoi(theString) < -10) {
							// Set our tag and sender
							//std::cout << "found a tag" << std::endl;
							currTag = lexical_cast<int>(*v);
							v+=2;
							currSender = lexical_cast<int>(*v);
							v-=2;

							//std::cout << "Tag: " << lexical_cast<std::string>(*v) << std::endl;

							if ( (tag == -1 || currTag == tag) && (senderID == -1 || senderID == currSender) ) {
								numMessages += 1;
								//std::cout << "incrementing numMessages" << std::endl;
							}
							else {
								//v += (s->second.size());
								break;
							}
							// If we find a tag (any number less than 10) we need to create a new vector
							// element within the message holder vector
							messageHolder.push_back(std::vector<std::string>());
						}
						//std::cout << messageHolder.size() << " Adding message: " << lexical_cast<std::string>(*v) << std::endl;
						// Add the item to the message holder
						messageHolder[messageHolder.size() - 1].push_back(lexical_cast<std::string>(*v));

						// Erase the element if eraseMessage is true. Otherwise, increment
						// the iterator
						if (eraseMessage) {
							//std::cout << "erasing " << *v << std::endl;
							v = s->second.erase(v);
						}
						else {
							++v;
						}
					}
					// If are array is empty, we should erase it
					/*if (s->second.size() == 0) {
						//std::cout << "erasing map pair" << std::endl;
						s = swarmMap_->erase(s);
					}
					else {
						++s;
					}*/
					++s;
				}
				else {
					++s;
				}
			}
			if (block && numMessages < 1) {
				usleep(2500);
			}
			else {
				//std::cout << "going to break" << std::endl;
				break;
			}
		}
		//std::cout << "out of loop?" << std::endl;
	}
	return numMessages;
}

void putArrayInSHM(std::vector<std::string> theArray) {

}
