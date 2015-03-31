
import Sire.Stream
import Sire.Qt

t = Sire.Qt.QTime()

filename = "test/SireStream/tmp_testdata.sire"

t.start()
header = Sire.Stream.getDataHeader(filename)

print(("Getting the header took %d ms" % t.elapsed()))

print((header.dataType()))
print((header.requiredLibraries()))

print((header.createdBy()))
print((header.createdWhere()))

print((header.requiredMemory()))
print((header.compressionRatio()))
print((header.digest()))
print((header.repository()))
print((header.buildVersion()))
print((header.systemInfo()))

print("Loading the system...")
t.start()
system = Sire.Stream.load(filename)

print((system.energies()))

print(("Loading the first time took %d ms" % t.elapsed()))

print("\nLoading the system again...")
t.start()
system = Sire.Stream.load(filename)

print(("Loading the second time took %d ms" % t.elapsed()))
