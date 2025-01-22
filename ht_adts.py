class HashTableIncrement(dict):
    def add(self, key: str, value: int =1):
        incremented_value = self.get(key, 0) + value
        self.update({key:incremented_value})


class HashTableAppend(dict):
    def add(self, key: str, value: any):
        appended_value = self.get(key, [])
        appended_value.append(value)
        self.update({key:appended_value})


class HashTableSet(dict):
    def add(self, key: str, location: int|str):
        value_set = self.get(key, set())
        value_set.add(location)
        self.update({key:value_set})


class HashTableString(dict):
    def add(self, key: str, string_addition: str):
        self.update({key:self.get(key, '') + string_addition})

