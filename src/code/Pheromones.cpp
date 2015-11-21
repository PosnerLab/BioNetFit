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
	if (swarm_->options.useCluster) {
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
					//TODO: Destructrr gets called when we go out of scope. Let's call it somewhere else..
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
			mutex_ = segment_->construct<interprocess_mutex>("MyMutex")();

			vecString_ = mystring;
			swarmMap_ = mymap;
			swarmVec_ = myshmvec;
		}
		else {
			MyShmString *mystring = new MyShmString(charallocator);
			MyMap_ *mymap = segment_->find<MyMap_>("SwarmMap").first;
			MyVector_ *myshmvec = segment_->find<MyVector_>("myvec").first;
			mutex_ = segment_->find<interprocess_mutex>("MyMutex").first;

			vecString_ = mystring;
			swarmMap_ = mymap;
			swarmVec_ = myshmvec;
		}
	}
}

void Pheromones::sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message) {
	// Using MPI
	if (swarm_->options.useCluster) {

		std::vector<int> receivers;

		// Sending to the entire swarm
		if (receiverID == -1) {
			std::vector<std::string> runningParticles;

			// First we need to get a list of all running particles from the master
			world_->send(0, GET_RUNNING_PARTICLES, "");
			world_->recv(0, SEND_RUNNING_PARTICLES, runningParticles);

			//for (auto p: runningParticles) {
			for (auto p = runningParticles.begin(); p != runningParticles.end(); ++p) {
				receivers.push_back(stoi(*p));
			}
		}
		else {
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		if (tag == -1) {
			tag = mpi::any_tag;
		}

		// Loop through receivers and perform the send operation
		for (std::vector<int>::iterator i = receivers.begin(); i != receivers.end(); ++i) {
			// Blocking send
			if (block) {
				world_->send(receiverID, tag, message);
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

			{
				scoped_lock<interprocess_mutex> lock(*mutex_);

				*vecString_ = std::to_string(static_cast<long long int>(GET_RUNNING_PARTICLES)).c_str(); // Set tag
				swarmVec_->push_back(*vecString_); // Push tag
				*vecString_ = std::to_string(static_cast<long long int>(rand())).c_str(); // Set unique ID
				swarmVec_->push_back(*vecString_); // Store unique ID
				*vecString_ = std::to_string(static_cast<long long int>(senderID)).c_str(); // Set sender
				swarmVec_->push_back(*vecString_); // Push sender
				*vecString_ = std::to_string(static_cast<long long int>(0)).c_str(); // Set message size of 0
				swarmVec_->push_back(*vecString_); // Store message size

				swarmMap_->insert(std::pair<const int,MyVector_>(0,*swarmVec_)); // Insert vector to first position (master slot) of map
				swarmVec_->clear(); // Clear our message vector
			}

			// Receive a message from master containing list of running particles
			//std::vector<std::vector<std::string>> messageHolder;
			recvMessage(0,senderID,SEND_RUNNING_PARTICLES,true, univMessageReceiver);
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
			//for (std::vector<std::string>::iterator r = messageHolder[0].begin(); r != messageHolder[0].end(); ++r) {
			//	receivers.push_back( atoi(r->c_str()) );
			//}
			swarmMsgHolderIt sm = univMessageReceiver.find(SEND_RUNNING_PARTICLES);
			for (std::vector<std::string>::iterator m = sm->second.message.begin(); m != sm->second.message.end(); ++m) {
				receivers.push_back( atoi(m->c_str()) );
				//std::cout << "adding receiver: " << m->c_str() << std::endl;
			}
		}
		else {
			// If we're not sending to swarm, add receiverID to list of message receivers
			// In this case, the list will only contain 1 item
			//std::cout << "adding receiverr: " << receiverID << std::endl;
			receivers.push_back(receiverID);
		}

		unsigned long int messageID = rand();

		{
			scoped_lock<interprocess_mutex> lock(*mutex_);

			// Add tag, sender, message size, and message(s) if applicable
			for(std::vector<int>::iterator r = receivers.begin(); r != receivers.end(); ++r) {
				int numMsgElements = message.size();

				*vecString_ = std::to_string(static_cast<long long int>(tag)).c_str(); // Set tag
				swarmVec_->push_back(*vecString_); // Store tag

				*vecString_ = std::to_string(static_cast<long long int>(messageID)).c_str(); // Set unique ID
				swarmVec_->push_back(*vecString_); // Store unique ID

				*vecString_ = std::to_string(static_cast<long long int>(senderID)).c_str(); // Set sender
				swarmVec_->push_back(*vecString_); // Store sender

				*vecString_ = std::to_string(static_cast<long long int>(numMsgElements)).c_str(); // Set message size
				swarmVec_->push_back(*vecString_); // Store message size

				// Insert messages in vector
				for(std::vector<std::string>::iterator m = message.begin(); m != message.end(); ++m) {
					*vecString_ = m->c_str();
					swarmVec_->push_back(*vecString_);
					//std::cout << "inserting " << *m << std::endl;
				}

				//std::cout << "storing vector in swarm map to receiver: " << *r << " with tag: " << tag << std::endl;
				// Insert vector to map
				swarmMap_->insert(std::pair<const int,MyVector_>(*r,*swarmVec_));
				//string tests = lexical_cast<string>(*(swarmMap_->find(*r)->second.begin()));
				//std::cout << "swarmmap is now: " << swarmMap_->size() << " " <<  tests << std::endl;
				swarmVec_->clear();
			}
		}

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
}

//int Pheromones::recvMessage(signed int senderID,const int receiverID, int tag, bool block, std::vector<std::vector<std::string>> &messageHolder, bool eraseMessage) {
int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage) {
	int numMessages = 0;
	if (swarm_->options.useCluster) {
		std::vector<std::string> mpiMessage;

		if (senderID == -1) {
			senderID = mpi::any_source;
		}

		if (tag == -1) {
			tag = mpi::any_tag;
		}

		while (1) {
			// TODO: Need to put this in a loop to gather more than one message at a taime
			if (block) {
				recvStatus_ = world_->recv(senderID, tag, mpiMessage);
			}
			else {
				recvRequest_ = world_->irecv(senderID, tag, mpiMessage);
			}

			// If we have any messages in our messageHolder, let's proces them
			if (messageHolder.size()) {

				// If we are not blocking (used irecv), need to wait on the message to receive its metadata
				if (!block) {
					recvStatus_ = recvRequest_.wait();
				}

				// Construct our swarmMessage and fill it with message and metadata
				swarmMessage smessage;
				//for (auto m: mpiMessage) {
				for (auto m = mpiMessage.begin(); m != mpiMessage.end(); ++m) {
					smessage.message.push_back(*m);
					smessage.tag = recvStatus_.tag();
					smessage.sender = recvStatus_.source();
				}

				// Insert the message into our message holder and increment numMessages
				messageHolder.insert(std::pair<int, swarmMessage>(tag, smessage));
				numMessages++;
			}

			// If we didn't receive any messages, time to clean up and exit the loop
			else if (!block) {
				recvRequest_.cancel();
				break;
			}
			else {
				break;
			}
		}
	}
	else {
		//std::cout << "in rcv msg" << std::endl;

		while (1) {
			if (!swarmMap_->empty()) {
				//std::cout << "sm size: " << swarmMap_->size() << std::endl;

				MyMap_::iterator s = swarmMap_->find(receiverID);

				if (s != swarmMap_->end()) {
					//std::cout << receiverID << " found message" << std::endl;
					// Lock the mutex so no other processes can access the shared memory
					scoped_lock<interprocess_mutex> lock(*mutex_);

					int currTag;
					int currSender;

					for (MyVector_::iterator v = s->second.begin(); v != s->second.end(); ) {
						std::string theString = lexical_cast<std::string>(*v);
						//std::cout << "msg: " << *v << std::endl;
						if (!theString.empty() && stoi(theString) < -10) {
							swarmMessage smessage;

							currTag = lexical_cast<int>(*v);
							v+=2;
							currSender = lexical_cast<int>(*v);

							if ( (tag == -1 || currTag == tag) && (senderID == -1 || senderID == currSender) ) {
								numMessages += 1;
								smessage.tag = currTag;
								//std::cout << "currTag is " << currTag << std::endl;
								smessage.sender = currSender;
								//std::cout << "currSender is " << currSender << std::endl;
								--v;
								smessage.id = lexical_cast<int>(*v);
								//std::cout << "id " << smessage.id << std::endl;

								v+=2;
								int size = lexical_cast<int>(*v);
								//std::cout << "size is " << size << std::endl;

								for (int i = 0; i < size; ++i) {
									++v;
									smessage.message.push_back(lexical_cast<std::string>(*v));
								}
								//std::cout << "inserting: " << lexical_cast<std::string>(*v) << std::endl;
								messageHolder.insert(std::pair<int, swarmMessage>(currTag, smessage));
								//messageHolder[currTag] = smessage;

								if (eraseMessage) {
									// Jump to beginning of message
									v -= (size+3);
									for (int i = 0; i <= (size+3); ++i) {
										//std::cout << "erasing: " << *v << std::endl;
										//sleep(1);
										v = s->second.erase(v);
									}
									//sleep(10);
								}
								else {
									++v;
								}
							}
						}
					}
					if (s->second.empty()) {
						//std::cout << "erasing map pair" << std::endl;
						s = swarmMap_->erase(s);
					}
				}
			}
			if (block && numMessages < 1) {
				usleep(10000);
				//std::cout << "wait" << std::endl;
			}
			else {
				//std::cout << "going to break" << std::endl;
				break;
			}
		}
	}

	/*
		while (1) {
			if (!swarmMap_->empty()) {
				//std::cout << "sm size: " << swarmMap_->size() << std::endl;

				MyMap_::iterator s = swarmMap_->find(receiverID);

				if (s != swarmMap_->end()) {
					scoped_lock<interprocess_mutex> lock(*mutex_);

					//std::cout << "first: " << s->first << " second: " << *(s->second.begin()) << std::endl;
					int currTag;
					int currSender;
					// TODO: Currently we can't have negative numbers as part of a message, since they will be confused as tags.  This needs fixed.
					for (MyVector_::iterator v = s->second.begin(); v != s->second.end(); ) {
						//std::cout << receiverID << ": " << *v << std::endl;
						std::string theString = lexical_cast<std::string>(*v);
						// TODO: We have to lexical cast to string, THEN convert to int because a lexical cast to int will throw an exception if the input can't be easily converted to an integer.  Is there a more elegant way to do this?
						if (!theString.empty() && stoi(theString) < -10) {
							// Set our tag and sender

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
					// If are array is empty, we should remove the int/array pair from the map
					if (s->second.empty()) {
						//std::cout << "erasing map pair" << std::endl;
						s = swarmMap_->erase(s);
					}
					//else {
					//	++s;
					//}
					//++s;
				}
				//else {
				//	++s;
				//}
			}
			if (block && numMessages < 1) {
				usleep(10000);
				//std::cout << "wait" << std::endl;
			}
			else {
				//std::cout << "going to break" << std::endl;
				break;
			}
		}
	 */
	return numMessages;
}
