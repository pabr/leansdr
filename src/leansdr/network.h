// This file is part of LeanSDR Copyright (C) 2016-2022 <pabr@pabr.org>.
// See the toplevel README for more information.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef LEANSDR_NETWORK_H
#define LEANSDR_NETWORK_H

#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

namespace leansdr {

// [udp_output] sends over UDP, one item per packet.

template<typename T>
struct udp_output : runnable {
  udp_output(scheduler *sch, pipebuf<T> &_in, const char *udpaddr,
	     int _batch_size)
    : runnable(sch, _in.name),
      in(_in),
      batch_size(_batch_size)
  {
    sock = socket(AF_INET, SOCK_DGRAM, 0);
    if ( sock < 0 ) fatal("socket");
    struct sockaddr_in addr;
    memset(&addr, 0, sizeof(addr));
    addr.sin_family = AF_INET;
    addr.sin_port = ntohs(0);
    addr.sin_addr.s_addr = INADDR_ANY;
    if ( bind(sock,(sockaddr*)&addr,sizeof(addr)) < 0 ) fatal("bind");
    addr.sin_family = AF_INET;
    {
      const char *sep = strchr(udpaddr, ':');
      if ( ! sep ) fail("Expected IP:PORT");
      int port = atoi(sep+1);
      addr.sin_port = ntohs(port);
      char *ipaddr = strndup(udpaddr, sep-udpaddr);
      if ( sch->verbose )
	fprintf(stderr, "Sending UDP to %s:%d\n", ipaddr, port);
      int res = inet_aton(ipaddr, &addr.sin_addr);
      free(ipaddr);
      if ( ! res ) fatal("inet_aton");
    }
    if ( connect(sock,(sockaddr*)&addr,sizeof(addr)) < 0 ) fatal("connect");
  }
  void run() {
    while ( in.readable() >= batch_size ) {
      // Return values for asynchronous protocols are not very useful.
      // Just ignore them.
      (void)send(sock, in.rd(), sizeof(T)*batch_size, 0);
      in.read(batch_size);
    }
  }
  ~udp_output() {
    close(sock);
  }
 private:
  pipereader<T> in;
  int sock;
  int batch_size;
};

}  // namespace

#endif  // LEANSDR_NETWORK_H
