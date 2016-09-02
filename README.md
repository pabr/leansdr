leansdr: Lightweight, portable software-defined radio.
Copyright (C) 2016 <pabr@pabr.org>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


**leansdr** consists of:
* A simple data-flow framework for signal processing
* A library of software-defined radio functions
* Applications built on top of the above.

Currently the main application is **leandvb**.

# leandvb

**leandvb** is a DVB-S demodulator designed for speed rather
than sensitivity.  See http://www.pabr.org/radio/leandvb .

## Quick start guide

```
git clone http://github.com/pabr/leansdr.git
cd leansdr/src/apps
make
```

### Receiving DATV transmissions from the ISS with a RTL-SDR dongle:

```
rtl_sdr  -f $DOWNCONVERTED_FREQ  -s 2400000  capture.iq
./leandvb  -f 2400e3  --sr 2000e3  --cr 1/2   < /tmp/capture.iq  > /tmp/capture.ts
mplayer capture.ts
```

### Troubleshooting

```
./leandvb_gui  --gui  -v  -d  -f 2400e3  --sr 2000e3  --cr 1/2  < /tmp/capture.iq  > /tmp/capture.ts
```

#### Live receiver with auto-detection of symbol rate and code rate:

```
rtl_sdr  -f $DOWNCONVERTED_FREQ  -s 2400000  -  |  ./leansdrscan  -v  ./leandvb_gui --gui  -f 2400e3  --sr 2000e3,1000e3,500e3,250e3  --cr 1/2,2/3,3/4,5/6,7/8  -  |  mplayer  -cache 128  -
```
