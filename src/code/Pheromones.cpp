/*============================================================================
// Name        : Pheromones.cpp
// Authors     : Brandon Thomas, Abolfazl Razi
// Version     : 2.0
// Lat Update: : 2017-01-15
// Copyright   :
// Description :
//============================================================================*/
#include "Pheromones.hh"
#include "Swarm.hh"

namespace mpi = boost::mpi;
using namespace boost::interprocess;
using boost::lexical_cast;


// TODO: Need lots of error checking here...

Pheromones::Pheromones() {
	swarm_ = 0;
	env_ = 0;
	world_ = 0;
}

Pheromones::~Pheromones() {
	if (swarm_->options.useCluster) {
	    int world_rank;
	    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	    if (swarm_->options.verbosity>=3) {
	    	cout << "my rank is " << world_rank << " starting ~environment()" << endl;
	    }
		//env_->~environment();
		MPI_Finalize();
	    if (swarm_->options.verbosity>=3) {

	    	cout << "~environment() done" << endl;
	    }

	}
	else {
		for (auto mq = smq_.begin(); mq != smq_.end(); ++mq) {
			(*mq)->~message_queue();
		}
	}
}

void Pheromones::init(Swarm *s) {
	swarm_ = s;
	unsigned int nSubPar= swarm_->options.swarmSize * swarm_->options.models.size()+1; //razi:number of subparticles, +1 since starts from 1

	// Using MPI
	if (swarm_->options.useCluster && (swarm_->options.clusterSoftware == "mpi" || swarm_->options.clusterSoftware == "slurm") ) {
		// Set up our MPI environment and communicator

		//std::cout<<"Pheromones Initialization: MPI\n";
		//env_ = new mpi::environment();
	    // Initialize the MPI environment
	    MPI_Init(NULL, NULL);
	    // Get the number of processes
	    int world_size;
	    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	    if (swarm_->options.verbosity>=3) {

	    	cout << "Defined mpi environment" << endl;
	    }
//		world_ = new mpi::communicator();
		 // Get the rank of the process
		    int world_rank;
		    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

		    //cout << "My rank is " << world_rank << " and I have just started." << endl;
		    if (swarm_->options.verbosity>=3) {

		    	cout << "Defined world communicator" << endl;
		    }
		//  boost::mpi::environment env;
		 // boost::mpi::communicator world_;
		//scout << "Hello World! from process " << world_.rank() << endl;
	}
	else if (swarm_->options.clusterSoftware == "BNF2mpi"){
	    if (swarm_->options.verbosity>=3) {

			cout << "Detected BNF2mpi in Pheromones init()" << endl;
	    }
			//env_ = new mpi::environment();
		    MPI_Init(NULL, NULL);
			    // Get the number of processes
			    int world_size;
			    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
			    if (swarm_->options.verbosity>=3) {

			    	cout << "Defined mpi environment" << endl;
			    }
		//		world_ = new mpi::communicator();
				 // Get the rank of the process
				    int world_rank;
				    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

				    if (swarm_->options.verbosity>=3) {

				    	cout << "My rank is " << world_rank << "and I have just started." << endl;

				    	cout << "Defined world communicator" << endl;
				    }
				//  boost::mpi::environment env;
	}
	// Using IPC
	else {

		//std::cout << "init ipc" << std::endl;
		if (s->isMaster) {
			std::cout<<"Pheromones Initialization: IPC: Master\n";
			for (unsigned int i = 0; i <= nSubPar ; ++i) {
				//std::cout << "creating: " << std::toString(static_cast<long long int>(i)) << std::endl;
				message_queue *smq;

				try {
					smq = new message_queue(create_only, toString(i).c_str(), 100, 1000);
				} catch (boost::interprocess::interprocess_exception e) {
					message_queue::remove(toString(i).c_str());
					smq = new message_queue(create_only, toString(i).c_str(), 100, 1000);
				}
				smq_.push_back(smq);
			}
		}
		else {
			std::cout<<"Pheromones Initialization: IPC: Slave\n";
			for (unsigned int i = 0; i <= nSubPar; ++i) {
				message_queue *smq = new message_queue(open_only, toString(i).c_str());
				smq_.push_back(smq);
				//std::cout << "opening: " << std::toString(static_cast<long long int>(i)) << " with max size of " << smq_[i]->get_max_msg() << std::endl;
			}
		}
	}

	swarm_->commInit = true;
}

void Pheromones::sendToSwarm(int senderID, signed int receiverID, int tag, bool block, std::vector<std::string> &message, int messageID) {
	// Using MPI
	if (swarm_->options.useCluster) {
		std::vector<int> receivers;

		// Find out rank, size
		int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
		int world_size;
		MPI_Comm_size(MPI_COMM_WORLD, &world_size);


		// Construct the swarmMessage
		swarmMessage smessage;
		smessage.sender = senderID;


		// Sending to the entire swarm
		if (receiverID == -1) {
			//std::cout << "sending to entire swarm because receiver is -1" << std::endl;

			// TODO: Replace this exchange with a world_->sendrecv()
			smessage.tag = toString(GET_RUNNING_PARTICLES);
		    if (swarm_->options.verbosity>=3) {

		    	std::cout << "trying to get list of running particles.." << std::endl;
			// First we need to get a list of all running particles from the master
		    }

			std::string serializedMessage = serializeSwarmMessage(smessage);
		    if (swarm_->options.verbosity>=3) {

		    	cout << "RRR sending serializedMessage GET RUNNING PARTICLES 11 = " << serializedMessage << endl;
		    }
			//world_->send(0, GET_RUNNING_PARTICLES, serializedMessage);
			MPI_Send(serializedMessage.c_str(), serializedMessage.length(), MPI_CHAR, 0, GET_RUNNING_PARTICLES, MPI_COMM_WORLD);

			//std::cout << "trying to receive list of running particles.." << std::endl;
			//usleep(10000);

			// This swarmMessage will hold the list of running particles received from master
			serializedMessage.clear();
			//world_->recv(0, SEND_RUNNING_PARTICLES, serializedMessage);
			char buffer[256];
			MPI_Status status;

			MPI_Recv(buffer, serializedMessage.length(), MPI_CHAR, 0, SEND_RUNNING_PARTICLES, MPI_COMM_WORLD, &status);
		    if (swarm_->options.verbosity>=3) {

		    	std::cout << "received list of running particles.." << std::endl;
		    }
			swarmMessage rsmessage = deserializeSwarmMessage(std::string(buffer));

			//for (auto p: runningParticles) {
			for (auto p = rsmessage.message.begin(); p != rsmessage.message.end(); ++p) {
			    if (swarm_->options.verbosity>=3) {

			    	std::cout << "adding receiver: " << *p << std::endl;
			    }
				receivers.push_back(stoi(*p));
			}

		}
		else {
		    if (swarm_->options.verbosity>=3) {

		    	cout << "RRR sendToSwarm senderID: " << senderID << " receiverID: " << receiverID << " tag: " << tag << endl;

				std::cout << "Sending tag " << tag << " to: " << receiverID << std::endl;
		    }
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		smessage.tag = toString(tag);

		//for (auto m: mpiMessage) {
		for (auto m = message.begin(); m != message.end(); ++m) {
			//std::cout << "adding: " << *m << std::endl;
			smessage.message.push_back(*m);
		}

		std::string smString;
		{
		smString = serializeSwarmMessage(smessage);
		}
		const char *serializedMessage = smString.c_str();
	    if (swarm_->options.verbosity>=3) {

	    	cout << "send message is serializedMessage = " << serializedMessage << " or smString = " << smString << endl;
	    }
		// Loop through receivers and perform the send operation
		for (std::vector<int>::iterator i = receivers.begin(); i != receivers.end(); ++i) {
			// Blocking send
		    int world_size;
		    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

			if(receiverID < world_size){

				//usleep(10000);


				if (block) {
				    if (swarm_->options.verbosity>=3) {

				    	std::cout << "attempting a block send from " << senderID << " to " << receiverID << std::endl;
				    }
					try{
//					world_->send(receiverID, convertedTag, serializedMessage, smString.length());
						//world_->send(receiverID, convertedTag, smString);
						MPI_Send(smString.c_str(), smString.length(), MPI_CHAR, receiverID, tag, MPI_COMM_WORLD);




					} catch (interprocess_exception& e) {
						cout << e.what( ) << std::endl;
			    	}
				    if (swarm_->options.verbosity>=3) {

				    	std::cout << "block send from " << senderID << " to " << receiverID << " succeeded" << std::endl;
				    }
				}
				// Non-blocking send
				else {
				    if (swarm_->options.verbosity>=3) {

				    	std::cout << "attempting a non-block send from " << senderID << " to " << receiverID << std::endl;
				    }
					//MPI_Request req;
					//recvRequest_ = world_->isend(receiverID, convertedTag, serializedMessage, smString.length());
					//recvRequest_ = world_->isend(receiverID, convertedTag, smString);
					//MPI_Isend(smString.c_str(), smString.length(), MPI_CHAR, receiverID, tag, MPI_COMM_WORLD, &req);

					   //for (int p = smString.find("\n"); p != (int) std::string::npos; p = smString.find("\n"))
						//   smString.erase(p,1);


					int result = MPI_Send(smString.c_str(), smString.length(), MPI_CHAR, receiverID, tag, MPI_COMM_WORLD);
	                if (result == MPI_SUCCESS){
	            	    if (swarm_->options.verbosity>=3) {

	                	std::cout << "Rank " <<  world_rank << " sent message to " << receiverID << " tag: " << tag << " message: " << smString.c_str() << std::endl;
	            	    }
	            	}else{
	            	    if (swarm_->options.verbosity>=3) {

	            	    	cout << "ERROR MESSAGE NOT SENT!!!" << endl;
	            	    }
					}
					//MPI_Send(smString.c_str(), smString.length(), MPI_CHAR, receiverID, tag, MPI_COMM_WORLD);
	        	    if (swarm_->options.verbosity>=3) {

	        	    	std::cout << "non-block send from " << senderID << " to " << receiverID << " succeeded" << std::endl;
	        	    }
					//recvStatus_ = recvRequest_.wait();
					//std::cout << "tag: " << recvStatus_.tag() << std::endl;
					//std::cout << "error: " << recvStatus_.error() << std::endl;
				}

			}
		}
	}else {
		// Using IPC
		std::vector<int> receivers;

		// Construct the swarmMessage
		swarmMessage smessage;
		smessage.sender = senderID;

		int nPar = 0; //Raquel: added a counter for debugging

		// Sending to the entire swarm
		if (receiverID == -1) {
			//std::cout << "sending to entire swarm because receiver is -1" << std::endl;
			// First we need to get a list of all running particles from the master

			// Set tag to tell master we need a list of running particles
			smessage.tag = toString(GET_RUNNING_PARTICLES);
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
				nPar++; //Raquel: added a counter for debugging

			}

			cout << "RAQUEL: Found " << nPar << " particles running." << endl; //Raquel: added for debugging

		}
		else {
			//std::cout << "Sending to: " << receiverID << std::endl;
			// If we're not sending to the entire swarm, put only the target pID into the receivers list
			receivers.push_back(receiverID);
		}

		// Set the message tag as specified by the sender
		smessage.tag = toString(tag);

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
			std::cout << "### " << senderID << " sending " << smessage.tag << " to " << *i	<< ". ser: " << serializedMessage.data() << std::endl;
			//smq_[*i]->send(serializedMessage.data(), serializedMessage.size(), 0);
			//Raquel: added for debbuging
			bool sendAttempt2 = smq_[*i]->try_send(serializedMessage.data(), serializedMessage.size(), 0);


			if (sendAttempt2 == false){
				cout << "RAQUEL: queue was full when message was sent" << endl;

			}else{ cout << "RAQUEL: sending success" << endl; }
			
			//std::cout << "sent" << std::endl;
			if (block) {
				bool foundMessage = true;
				while (foundMessage) {
					//usleep(150000);
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

int Pheromones::recvMessage(signed int senderID, const int receiverID, int tag, bool block, swarmMsgHolder &messageHolder, bool eraseMessage, int messageID) {
	int numMessages = 0;
	std::string serializedMessage; //Raquel: moved this from inside the loop to here so it's not redeclared
	serializedMessage.resize(1000);
	swarmMessage smessage;

	if (swarm_->options.useCluster) {
		//swarmMessage smessage;
		//std::string serializedMessage;
		//cout << "RRR recvMessage senderID: " << senderID << " receiverID: " << receiverID << " tag: " << tag <<endl;
		std::string msg;
		MPI_Status status;
		int flag = 0;
		//char buffer[256];
		serializedMessage.resize(1000);

		while (1) {
			//std::cout << "rcv loop" << std::endl;
			//usleep(10000);

			MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
			//MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			//if(status == MPI_SUCCESS){
			if (flag==1) {
				//cout << "status tag: " << recvStatus->tag() << std::endl;
				//cout << "status source: " << recvStatus->source() << std::endl;
				//boost::optional<int> msgLength = recvStatus->count<char>();
				//char smChar[*msgLength];
				//char smChar[1000];
				//cout << "recvStatus->tag() " << recvStatus->tag() << " DONE_BREEDING " << DONE_BREEDING << endl;
				int count;

				MPI_Get_count(&status, MPI_CHAR, &count);


				int* pVal = &count;
				int l = *pVal;

				char *buffer = new char[l];
			    if (swarm_->options.verbosity>=3) {

			    	cout << "message size " << l << endl;
			    }
				//if (block) {
			    if (swarm_->options.verbosity>=3) {

			    	cout << "trying a blocking receive to receiverID " << receiverID << " from senderID " << senderID  << endl;
			    }
				//cout << "Tag: " << tag << "; smChar length: " << sizeof(smChar) << "; *msgLength: " << *msgLength << endl;
				//world_->recv(senderID, tag, smChar, *msgLength);
				//world_->recv(senderID, recvStatus->tag(), smChar, *msgLength);
					//world_->recv(senderID, recvStatus->tag(), msg);
				MPI_Recv(buffer, l, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
				msg = std::string(buffer);
				//cout <<"RRR received message is smChar = " << smChar << endl;
			    if (swarm_->options.verbosity>=3) {

			    	cout <<"RRR received buffer is = " << buffer << endl;
			    	cout <<"RRR received message is = " << msg << endl;

					cout << "Message size is " << msg.length() << endl;
			    }
				//world_->recv(boost::mpi::any_source, tag, *msgLength);
				block = false;
			    if (swarm_->options.verbosity>=3) {

			    	cout << "RRR passed receive step." << endl;
			    }

				//smessage = deserializeSwarmMessage(std::string(smChar));
				if(msg.length()>0){
					smessage = deserializeSwarmMessage(msg);


				}else{
				    if (swarm_->options.verbosity>=3) {

				    	cout << "MESSAGE EMPTY!!!" << endl;
				    }
					//smessage.tag = std::to_string(recvStatus->tag());
					//smessage.sender = senderID;
					//smessage.id = messageID;
				}
			    if (swarm_->options.verbosity>=3) {

			    	cout << "RRR deserialized message" << endl;
			    }
				msg.clear();
			    if (swarm_->options.verbosity>=3) {

			    	cout << "Cleaned serialized message" << endl;
			    }
				// If we have any messages in our messageHolder, let's process them
				//std::cout << "messageholder not empty: " << smessage.tag << ":" << smessage.sender << std::endl;

				// Make sure our message matches the sender, tag, and id we requested
				if ( (tag == -1 || tag == stoi(smessage.tag)) && (senderID == -1 || senderID == smessage.sender) && (messageID == -1 || messageID == smessage.id)) {
					// Insert the message into our message holder and increment numMessages
				    if (swarm_->options.verbosity>=3) {

				    	std::cout << "inserting message in the holder" << std::endl;
				    }
					messageHolder.insert(std::pair<int, swarmMessage>(stoi(smessage.tag), smessage));

					++numMessages;

					// If user doesn't want to erase the message, put it back in the queue
					if (!eraseMessage) {
					    if (swarm_->options.verbosity>=3) {

					    	std::cout << "Not erasing message" << std::endl;
					    }
						sendToSwarm(senderID, receiverID, stoi(smessage.tag), false, smessage.message);
					}

					// Clear out the smessage for next use
					clearSwarmMessage(smessage);
				}
				else {
				    if (swarm_->options.verbosity>=3) {

				    	std::cout << "putting it back in the queue..." << std::endl;
				    	std::cout << "@@@@@@@@@@@@@@@@@@@@@@@ Expected: " << tag << "; Received: " << smessage.tag << "@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
				    }
					sendToSwarm(smessage.sender, receiverID, stoi(smessage.tag), false, smessage.message);

					if(tag == DONE_BREEDING && stoi(smessage.tag) == SIMULATION_END ){

						if(swarm_->options.verbosity >= 3){
							std::cout << "SENDER REINSERTED = " << smessage.sender << endl;
						}
						swarm_->fixRunningParticle(smessage.sender);

						std::vector<unsigned int> NewlyFinishedParticles = swarm_->checkMasterMessages();

						for(auto i = NewlyFinishedParticles.begin(); i!=NewlyFinishedParticles.end(); ++i){
							if (swarm_->options.verbosity>=3) {
						    	std::cout << "@@@@@@@@@@@@@@@@@@@@@@@ PROCESSING LATE PAR: " << *i << std::endl;
						    	std::cout << "from pheromones: " << swarm_->subparticleCurrParamSets_.size() << std::endl;
						    }
							std::map<std::string,double> mysimParams;
							int mid = fcalcMID(*i, swarm_->options.models.size());
							int pID = fcalcParID(*i, swarm_->options.models.size());

							auto fp = swarm_->options.models.at(mid)->freeParams_.begin();

							for (unsigned int j = 0; j < swarm_->subparticleCurrParamSets_[pID][mid].size(); j++) {
//								std::map<unsigned int, std::map<unsigned int, std::vector<double> > > subparticleCurrParamSets_;   //razi: PID, MID, PARAMA VALUESS

								//paramSet.insert(pair<string, double> (fp->first, subparticleCurrParamSets_[bestP][mid][i])); //make this work later, to clone the best particle
								mysimParams.insert(pair<std::string, double> (fp->first, swarm_->subparticleCurrParamSets_[pID][mid][j]));

								++fp;
							}
							swarm_->processLateParticles(mysimParams, *i, true, 0);

						}



					}
				}

				break;
			}
			else if (!block) {
				break;
			}
		}
	}
	else {
		//cout << "RAQUEL entering recvMessage entering while loop" << endl;
		//cout << "RAQUEL receiver ID: " << receiverID << " senderID: " << senderID << endl;
		unsigned int priority = 0;

		while(1) {
			//cout << "RAQUEL inside recvMessage entered while loop" << endl;

			message_queue::size_type recvd_size;

			//std::cout << "rcv loop" << std::endl;
			if (block) {
//std::cout << "trying a blocking receive to " << receiverID << " from " << senderID << std::endl;

				smq_[receiverID]->receive(&serializedMessage[0], 1000, recvd_size, priority);


				block = false;
//std::cout << "blocking receive to " << receiverID << " from " << senderID << " succeeded with recv_size of " << recvd_size << " and message of " << serializedMessage << std::endl;
			}
			else {

//std::cout << "trying a non-blocking receive to " << receiverID << " from " << senderID << std::endl;
				bool hasMessage = smq_[receiverID]->try_receive(&serializedMessage[0], 1000, recvd_size, priority);

//std::cout << "non-blocking receive to " << receiverID << " from " << senderID << " succeeded with recv_size of " << recvd_size << " and bool of " << hasMessage << std::endl;
				if (!hasMessage) {
//std::cout << "breaking" << std::endl;
					break;
				}
			}
			// If we have any messages in our messageHolder, let's process them

			// De-serialize the smessage
			serializedMessage.resize(recvd_size);
			swarmMessage smessage = deserializeSwarmMessage(serializedMessage);

			serializedMessage.clear();
			serializedMessage.resize(1000);

sleep(1);

			//razi: uncomment later
			//std::cout << "smessage not empty: " << smessage.tag << std::endl;
			//std::cout << "tag: " << tag << ":" << smessage.tag << std::endl;
			//std::cout << "sender: " << tag << ":" << smessage.sender << std::endl;
			//std::cout << "id: " << tag << ":" << smessage.id << std::endl;

			// Make sure our message matches the sender, tag, and id we requested
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

				if (!block) {
					break;
				}
			}
		}
	}
	//cout << "RAQUEL inside recvMessage exiting with num messages: " << numMessages << endl;

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
	{
	boost::archive::text_oarchive oa(oss);
	oa << sm;
	}
	std::string serializedMessage(oss.str());

	return serializedMessage;
}

Pheromones::swarmMessage Pheromones::deserializeSwarmMessage(std::string sm) {


	swarmMessage smessage;

	//char ch = sm.back();
	//cout << "Last character @" << ch << "@" << endl;
	if(sm.empty()){
	//if(ch=='\0'){
		//sm.back() = '\n';
    	sm="\n";

	    if (swarm_->options.verbosity>=3) {

			cout << "String empty" << endl;
	    }
	}else{

		//cout << "last character is not empty, erasing" << endl;
		//sm.erase(sm.size()-1);
		//ch = sm.back();
		//cout << "new last character @" << ch << "@"<< endl;
	}


	if (swarm_->options.verbosity>=3) {

		cout << "serializing message start SM: " << sm << endl;
	}

	{
	std::stringstream iss;
	iss.str(sm);
	boost::archive::text_iarchive ia(iss);
	if (swarm_->options.verbosity>=3) {
		cout << "done serializing message" << endl;
	}
	ia >> smessage;
	}
	return smessage;
}

int Pheromones::getRank() {
	//return world_->rank();
    int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	return world_rank;

}
