# Environment setup

The RS-MEVC pipeline calls external binaries (NetMHCpan 4.1b, NetMHCIIpan 4.3i)
that must be installed locally. The refactored scripts receive the binary paths
via CLI arguments, so they do not require any shell-wide setup. The templates
in this directory are provided only as a convenience.

## NetMHCpan 4.1b (MHC-I)

```bash
cp scripts/setup/netmhcpan_setup.sh.example ~/netmhcpan_setup.sh
# edit NETMHCPAN_HOME inside ~/netmhcpan_setup.sh
source ~/netmhcpan_setup.sh
```

## NetMHCIIpan 4.3i (MHC-II)

```
cp scripts/setup/netmhciipan_setup.sh.example ~/netmhciipan_setup.sh
# edit NETMHCIIPAN_HOME inside ~/netmhciipan_setup.sh
source ~/netmhciipan_setup.sh
```