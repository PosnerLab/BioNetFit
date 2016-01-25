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

Pheromones::Pheromones() {
	swarm_ = 0;
	env_ = 0;
	world_ = 0;
}

void Pheromones::init(Swarm *s) {
	swarm_ = s;

	// Using MPI
	if (swarm_->options.useCluster) {
		// Set up our MPI environment and communicator
		env_ = new mpi::environment();
		world_ = new mpi::communicator();
	}
	// Using IPC
	else {

		//std::cout << "init ipc" << std::endl;
		if (s->isMaster) {
			for (unsigned int i = 0; i <= swarm_->options.swarmSize; ++i) {
				// TODO: The type conversion here is horrible

				//std::cout << "creating: " << std::to_string(static_cast<long long int>(i)) << std::endl;
				message_queue::remove(std::to_string(static_cast<long long int>(i)).c_str());
				message_queue *smq = new message_queue(create_only, std::to_string(static_cast<long long int>(i)).c_str(), 100, 1000);
				smq_.push_back(smq);
			}
		}
		else {
			for (unsigned int i = 0; i <= swarm_->options.swarmSize; ++i) {
				message_queue *smq = new message_queue(open_only, std::to_string(static_cast<long long int>(i)).c_str());
				smq_.push_back(smq);
				//std::cout << "opening: " << std::to_string(static_cast<long long int>(i)) << " with max size of " << smq_[i]->get_max_msg() << std::endl;
			}
		}
	}
}

void Pheromones::sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message, int messageID) {
	// Using MPI
	if (swarm_->options.useCluster) {
		std::vector<int> receivers;

		// Construct the swarmMessage
		swarmMessage smessage;
		smessage.sender = senderID;

		// Sending to the entire swarm
		if (receiverID == -1) {
			//std::cout << "sending to entire swarm because receiver is -1" << std::endl;

			// TODO: Replace this exchange with a world_->sendrecv()
			smessage.tag = std::to_string(static_cast<long long int>(GET_RUNNING_PARTICLES));
			//std::cout << "trying to get list of running particles.." << std::endl;
			// First we need to get a list of all running particles from the master

			std::string serializedMessage = serializeSwarmMessage(smessage);
			world_->send(0, GET_RUNNING_PARTICLES, serializedMessage);
			//std::cout << "trying to receive list of running particles.." << std::endl;

			// This swarmMessage will hold the list of running particles received from master
			serializedMessage.clear();
			world_->recv(0, SEND_RUNNING_PARTICLES, serializedMessage);
			//std::cout << "received list of running particles.." << std::endl;
			swarmMessage rsmessage = deserializeSwarmMessage(serializedMessage);

			//for (auto p: runningParticles) {
			for (auto p = rsmessage.message.begin(); p != rsmessage.message.end(); ++p) {
				//std::cout << "adding receiver: " << *p << std::endl;
				receivers.push_back(stoi(*p));
			}
		}
		else {
			//std::cout << "Sending " << tag << " to: " << receiverID << std::endl;
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		smessage.tag = std::to_string(static_cast<long long int>(tag));

		//for (auto m: mpiMessage) {
		for (auto m = message.begin(); m != message.end(); ++m) {
			//std::cout << "adding: " << *m << std::endl;
			smessage.message.push_back(*m);
		}

		std::string smString = serializeSwarmMessage(smessage);
		const char *serializedMessage = smString.c_str();
		//std::cout << "send message is: " << serializedMessage << std::endl;

		// Loop through receivers and perform the send operation
		for (std::vector<int>::iterator i = receivers.begin(); i != receivers.end(); ++i) {
			// Blocking send
			if (block) {
				//std::cout << "attempting a block send from " << senderID << " to " << receiverID << std::endl;
				world_->send(receiverID, tag, serializedMessage, smString.length());
				//std::cout << "block send from " << senderID << " to " << receiverID << " succeeded" << std::endl;
			}
			// Non-blocking send
			else {
				//std::cout << "attempting a non-block send from " << senderID << " to " << receiverID << std::endl;
				recvRequest_ = world_->isend(receiverID, tag, serializedMessage, smString.length());
				//std::cout << "non-block send from " << senderID << " to " << receiverID << " succeeded" << std::endl;
				//recvStatus_ = recvRequest_.wait();
				//std::cout << "tag: " << recvStatus_.tag() << std::endl;
				//std::cout << "error: " << recvStatus_.error() << std::endl;
			}
		}
	}

	else {
		// Using IPC
		std::vector<int> receivers;

		// Construct the swarmMessage
		swarmMessage smessage;
		smessage.sender = senderID;

		// Sending to the entire swarm
		if (receiverID == -1) {
			//std::cout << "sending to entire swarm because receiver is -1" << std::endl;
			// First we need to get a list of all running particles from the master

			// Set tag to tell master we need a list of running particles
			smessage.tag = std::to_string(static_cast<long long int>(GET_RUNNING_PARTICLES));
			//std::cout << "trying to get list of running particles.." << std::endl;

			// Serialize the message
			std::string serializedMessage = serializeSwarmMessage(smessage);
			serializedMessage.resize(1000);

			smq_[0]->send(serializedMessage.data(), sizeof(serializedMessage), 0);
			//std::cout << "trying to receive list of running particles.." << std::endl;

			// This swarmMessage will hold the list of running particles received from master

			message_queue::size_type recvd_size;

			// Receive and de-serialize the message
			std::stringstream iss;
			serializedMessage.clear();
			serializedMessage.resize(1000);
			unsigned int priority;

			smq_[0]->receive(&serializedMessage[0], 1000, recvd_size, priority);

			serializedMessage.resize(recvd_size);
			swarmMessage rsmessage = deserializeSwarmMessage(serializedMessage);

			//std::cout << "received list of running particles.." << std::endl;

			//for (auto p: runningParticles) {
			for (auto p = rsmessage.message.begin(); p != rsmessage.message.end(); ++p) {
				//std::cout << "adding receiver: " << *p << std::endl;
				receivers.push_back(stoi(*p));
			}
		}
		else {
			//std::cout << "Sending to: " << receiverID << std::endl;
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		// Set the message tag as specified by the sender
		smessage.tag = std::to_string(static_cast<long long int>(tag));

		// Add the message array to the swarmMessage
		//for (auto m: mpiMessage) {
		for (auto m = message.begin(); m != message.end(); ++m) {
			//std::cout << senderID << " adding: " << *m << std::endl;
			smessage.message.push_back(*m);
		}

		// Set a random messageID
		if (messageID == -1) {
			messageID = rand();
		}

		smessage.id = messageID;

		// Serialize the swarmMessage
		std::string serializedMessage = serializeSwarmMessage(smessage);
		serializedMessage.resize(1000);

		// Loop through receivers and perform the send operation
		for (std::vector<int>::iterator i = receivers.begin(); i != receivers.end(); ++i) {
			//std::cout << senderID << " sending " << smessage.tag << " to " << *i	<< ". ser: " << serializedMessage.data() << std::endl;
			smq_[*i]->send(serializedMessage.data(), serializedMessage.size(), 0);
			//std::cout << "sent" << std::endl;
			if (block) {
				bool foundMessage = true;
				while (foundMessage) {
					usleep(150000);
					//std::cout << senderID << "loop" << std::endl;

					// TODO: This non-erasing recvMessage might slow down the receiver finding the message since it
					// requires a re-send every time we check

					// If we're blocking, we need to repeatedly check to see if the message still exists in the queue
					if (recvMessage(senderID, receiverID, tag, false, univMessageReceiver, false, messageID)) {
						foundMessage = true;
					}
					else {
						foundMessage = false;
					}
				}
			}
		}
	}
}

//int Pheromones::recvMessage(signed int senderID,const int receiverID, int tag, bool block, std::vector<std::vector<std::string>> &messageHolder, bool eraseMessage) {
int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {
	int numMessages = 0;
	if (swarm_->options.useCluster) {
		//swarmMessage smessage;
		std::string serializedMessage;

		while (1) {
			//std::cout << "rcv loop" << std::endl;

			if (boost::optional<boost::mpi::status> recvStatus = world_->iprobe(senderID, tag)) {
				//std::cout << "status: " << recvStatus->tag() << std::endl;

				boost::optional<int> msgLength = recvStatus->count<char>();
				char smChar[*msgLength+1];

				//if (block) {
				//std::cout << "trying a blocking receive to " << receiverID << " from " << senderID << std::endl;
				world_->recv(senderID, tag, smChar, *msgLength);
				block = false;

				swarmMessage smessage = deserializeSwarmMessage(std::string(smChar));
				serializedMessage.clear();

				// If we have any messages in our messageHolder, let's process them
				//std::cout << "messageholder not empty: " << smessage.tag << ":" << smessage.sender << std::endl;

				// Make sure our message matches the sender, tag, and id we requested
				if ( (tag == -1 || tag == stoi(smessage.tag)) && (senderID == -1 || senderID == smessage.sender) && (messageID == -1 || messageID == smessage.id)) {
					// Insert the message into our message holder and increment numMessages
					//std::cout << "inserting " << std::endl;
					messageHolder.insert(std::pair<int, swarmMessage>(stoi(smessage.tag), smessage));
					++numMessages;

					// If user doesn't want to erase the message, put it back in the queue
					if (!eraseMessage) {
						//std::cout << "Not erasing message" << std::endl;
						sendToSwarm(senderID, receiverID, tag, false, smessage.message);
					}

					// Clear out the smessage for next use
					clearSwarmMessage(smessage);
				}
				else {
					//std::cout << "putting it back in the queue..." << std::endl;
					sendToSwarm(senderID, receiverID, tag, false, smessage.message);
				}
			}
			else if (!block) {
				break;
			}
		}
	}
	else {
		// TODO: Move these to header so they don't get redeclared every time?
		std::string serializedMessage;
		serializedMessage.resize(1000);

		unsigned int priority = 0;

		while(1) {
			message_queue::size_type recvd_size;

			//std::cout << "rcv loop" << std::endl;
			if (block) {
				//std::cout << "trying a blocking receive to " << receiverID << " from " << senderID << std::endl;
				smq_[receiverID]->receive(&serializedMessage[0], 1000, recvd_size, priority);
				block = false;
				//std::cout << "blocking receive too " << receiverID << " from " << senderID << " succeeded with recv_size of " << recvd_size << " and message of " << serialized_string << std::endl;
			}
			else {
				//std::cout << "trying a non-blocking receive to " << receiverID << " from " << senderID << std::endl;
				bool hasMessage = smq_[receiverID]->try_receive(&serializedMessage[0], 1000, recvd_size, priority);

				//std::cout << "non-blocking receive to " << receiverID << " from " << senderID << " succeeded with recv_size of " << recvd_size << " and bool of " << hasMessage << std::endl;
				if (!hasMessage) {
					//std::cout << "breaking?" << std::endl;
					break;
				}
			}
			// If we have any messages in our messageHolder, let's process them
			//std::cout << "smessage not empty: " << serialized_string << std::endl;

			// De-serialize the smessage
			serializedMessage.resize(recvd_size);
			swarmMessage smessage = deserializeSwarmMessage(serializedMessage);


			serializedMessage.clear();
			serializedMessage.resize(1000);

			// Make sure our message matches the sender, tag, and id we requested
			/*
			std::cout << "tag: " << tag << ":" << smessage.tag << std::endl;
			std::cout << "sender: " << tag << ":" << smessage.sender << std::endl;
			std::cout << "id: " << tag << ":" << smessage.id << std::endl;
			 */
			if ( (tag == -1 || stoi(smessage.tag) == tag) && (senderID == -1 || senderID == smessage.sender) && (messageID == -1 || messageID == smessage.id)) {
				// Insert the message into our message holder and increment numMessages
				messageHolder.insert(std::pair<int, swarmMessage>(stoi(smessage.tag), smessage));
				//std::cout << "storing pair with tag of " << stoi(smessage.tag) << std::endl;
				++numMessages;

				// If user doesn't want to erase the message, put it back in the queue
				if (!eraseMessage) {
					//std::cout << "noerase" << std::endl;
					sendToSwarm(senderID, receiverID, tag, false, smessage.message);
				}
			}
			else {
				// If this isn't our message, put it back in the queue. This will go SLOW
				// unless we use a different queue for every particle

				//std::cout << "putting it back in the queue..." << std::endl;
				sendToSwarm(smessage.sender, receiverID, stoi(smessage.tag), false, smessage.message, smessage.id);
			}
		}
	}

	return numMessages;
}

void Pheromones::clearSwarmMessage(swarmMessage& sm) {
	sm.tag.clear();
	sm.id = 0;
	sm.sender = 0;
	sm.message.clear();
}

std::string Pheromones::serializeSwarmMessage(swarmMessage sm) {
	std::stringstream oss;
	boost::archive::text_oarchive oa(oss);
	oa << sm;
	std::string serializedMessage(oss.str());

	return serializedMessage;
}

Pheromones::swarmMessage Pheromones::deserializeSwarmMessage(std::string sm) {

	std::stringstream iss;
	iss.str(sm);
	boost::archive::text_iarchive ia(iss);
	swarmMessage smessage;
	ia >> smessage;

	return smessage;
}

int Pheromones::getRank() {
	return world_->rank();
}
