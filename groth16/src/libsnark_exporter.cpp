#include "libsnark_exporter.h"

Offset<Vector<Offset<KeyValue>>> make_configuration(FlatBufferBuilder &builder, vector<pair<string, string>> keyvalues)
{
  vector<Offset<KeyValue>> config;
  // add config key-value pairs
  for (auto kv = keyvalues.begin(); kv != keyvalues.end(); kv++)
  {
    config.emplace_back(CreateKeyValue(builder, builder.CreateString(kv->first),
                                       0, builder.CreateString(kv->second)));
  }
  return builder.CreateVector(config);
}