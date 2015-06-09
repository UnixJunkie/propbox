from __future__ import print_function, absolute_import

from .simple_futures import (new_future as _new_future,
                             new_future_exception as _new_future_exception)

# Propbox organizes dependencies between molecular property/descriptor
# calculations.

__version__ = "0.5"
__version_tuple = (0, 5, 0)

DEBUG = False

############### Propbox exceptions

class PropboxError(Exception):
    pass

# A PropboxKeyError is raised when the resolver does not implement a property
    
class PropboxKeyError(KeyError, PropboxError):
    """A PropboxKeyError is raised when the table does not implement a property name

    The two attributes are:
      `resolver` - the resolver under question
      `column_name` - the property name that does not exist

    """
    def __init__(self, resolver, column_name):
        super(PropboxKeyError, self).__init__(self, (resolver, column_name))
        self.resolver = resolver
        self.column_name = column_name
    def __repr__(self):
        return "PropboxKeyError(%r, %r)" % (self.resolver, self.column_name)
    def __str__(self):
        return str(self.column_name)
    
# A ResolverError describes why, for a row, it was not possible to get
# the given column name from the table.

class ResolverError(PropboxError):
    """A ResolverError is raised when a property cannot be computed.

    The three attributes are:
      `exception` - the exception describing the failure
      `table_name` - the name of the table where the failure occurred
      `column_name` - the property name that could not be resolved

    Note: the exception may be another ResolverError if the failure
    occurs further up the dependency chain. It may require several
    links to get to the fundamental exception.
    """
    
    def __init__(self, exception, table_name, column_name):
        self.exception = exception
        self.table_name = table_name
        self.column_name = column_name
    
    def get_error(self):
        """Internal helper function use to implement __str__

        Should not be used directly. (Let me know if you think
        it should be part of the public API.)
        """
        # Get the fully qualified name for the column
        if self.table_name is None:
            name = self.column_name
        else:
            name = "%s.%s" % (self.table_name, self.column_name)
        # Recursively expand to get a description of full
        # path that lead to the actual exception.
        if isinstance(self.exception, ResolverError):
            path, exception = self.exception.get_error()
            path = name + " -> " + path
        else:
            path = name
            exception = self.exception
        # return a path like "abc -> module.def -> mod2.ghi -> xyz"
        # and the exception that caused the error.
        return path, exception

    def get_original_exception(self):
        if isinstance(self.exception, ResolverError):
            return self.exception.get_original_exception()
        return self.exception
    
    def __str__(self):
        path, exception = self.get_error()
        return "%s: %r" % (path, exception)

##############################


### A "Resolver" adds 1 or more columns to a table.

class Resolver(object):
    resolver_name = None
    output_names = []  # List of column names that this knows about

    def resolve_column(self, name, table):
        # Add column 'column' to the table. May also add other columns.
        # Raises a PropboxKeyError if the name isn't supported
        raise NotImplementedError


# Some of the futures may contain an exception. Wrap them in a
# ResolveError exception so the error reporting has the right
# chain of descriptor names.
def wrap_future_exceptions(futures, table_name, column_name):
    new_futures = []
    for future in futures:
        exception = future.exception()
        if exception is not None:
            future = _new_future_exception(ResolverError(exception, table_name, column_name))
        new_futures.append(future)
    return new_futures



### A "Propbox" is a resolver which implements zero or more resolvers

class Propbox(Resolver):
    def __init__(self, resolvers=None, resolver_name=None):
        # All of the output descriptors from the resolvers
        self.output_names = []

        # Mapping from output name to the resolver which compute it.
        self._name_to_resolver = {}

        # Keep for later
        self._resolvers = []
        if resolvers is not None:
            for resolver in resolvers:
                self.add_resolver(resolver)

        self.resolver_name = resolver_name

    def __iter__(self):
        return iter(self._resolvers)

    def add_resolver(self, resolver):
        # Do one pass to ensure that the output names are unique
        # in the resolver, and unique in the propbox.
        column_names = set()
        for output_name in resolver.output_names:
            if output_name in column_names:
                raise ValueError(
                    "Resolver %r contains the output_name %r twice"
                    % (resolver, output_name))
            # Consider if resolver R1 implements A and B,
            # and resolver R2 implements B and C, and where
            # the value of B is different in R1 or R2.
            # Because a resolver can resolver more than one
            # column at a time, it's possible to get different
            # values for ["A", "B", "C"] than ["C", "B", "A"].
            #
            # This check prevents that problem from coming up.
            if output_name in self._name_to_resolver:
                raise ValueError(
                    "Resolver %r defines the output name %r, which was "
                    "already defined by resolver %r"
                    % (resolver, output_name,
                       self._name_to_resolver[output_name]))
            column_names.add(output_name)

        # Do another pass to update the registries
        for output_name in resolver.output_names:
            self._name_to_resolver[output_name] = resolver
            self.output_names.append(output_name)

        self._resolvers.append(resolver)

    def resolve_column(self, name, table):
        resolver = self._name_to_resolver.get(name, None)
        if resolver is None:
            raise PropboxKeyError(self, name)

        # Forward the request to the resolver
        resolver.resolve_column(name, table)

        if not table.has_computed(name):
            raise AssertionError(
                "Propbox %r resolver %r did not resolve the request for column name %r"
                % (self, resolver, name))

### An Alias adds a new way to refer to an exisiting name

class Aliases(Resolver):
    def __init__(self, alias_map, resolver_name=None):
        # output name -> input name
        self.alias_map = alias_map
        #self.input_names = sorted(alias_map.values())
        self.output_names = sorted(alias_map)
        self.resolver_name = resolver_name

    def resolve_column(self, name, table):
        try:
            alias_name = self.alias_map[name]
        except KeyError:
            raise PropboxKeyError(self, name)

        # Get the underlying name, wrap any exception with a description of
        # the alias step, and set the new results.
        futures = table.get_futures(alias_name)
        futures = wrap_future_exceptions(futures, table.table_name, name)
        table.set_futures(name, futures)

                
### A Module isolates a resolver into its own subtable


class ParentTableResolver(Resolver):
    def __init__(self, parent_table, resolver, input_map):
        self.parent_table = parent_table
        self.resolver = resolver
        self.input_map = input_map

        self.output_names = list(input_map) + resolver.output_names

    def resolve_column(self, name, subtable):
        if name in self.input_map:
            # Forward to the parent table
            parent_name = self.input_map[name]
            parent_futures = self.parent_table.get_futures(parent_name)

            # Wrap any exceptions with the alias information
            parent_futures = wrap_future_exceptions(parent_futures, table.table_name, name)

            # Copy the parent table values into the subtable
            table.set_futures(name, parent_futures)

        else:
            # Let the resolver handle it.
            self.resolver.resolve_column(name, table)


class Module(Resolver):
    def __init__(self, module_name, resolver, input_map, output_map):
        self.module_name = module_name # Used to fully quality a property name
        self.resolver = resolver
        
        # Map from subtable name to parent name. For example,
        #  {"my_mol": "mol"}
        # will map the subtable's request for "my_mol" into a
        # parent table request for "mol"
        self.input_map = input_map

        for name in resolver.output_names:
            if name in self.input_map:
                raise ValueError(
                    "The name %r exists in both the input_map (forwarded to %r) and in the resolver %r"
                    % (name, input_map[name], self.resolver))
        
        # Map from parent table name to subtable name. For example,
        #   {"BBB_v1": "prediction"}
        # will map the parent's request for "BBB_v1" into the
        # subtable's request for "prediction".
        self.output_map = output_map

        self.output_names = list(output_map)

        self.resolver_name = module_name

    def resolve_column(self, name, table):
        if name not in self.output_map:
            raise PropboxKeyError(self, name)
        
        subtable = table.get_cache_value(self.module_name)
        if subtable is None:
            # Figure out the subtable name
            if table.table_name is None:
                subtable_name = self.module_name
            else:
                subtable_name = table.table_name + "." + self.module_name

            # Get only the configuration items for this module
            config = {}
            prefix = self.module_name + "."
            n = len(prefix)
            for k, v in table.config.iteritems():
                if k[:n] == prefix:
                    config[k[n:]] = v
                    
            subtable = PropertyTable(
                resolver = ParentTableResolver(table, self.resolver, self.input_map),
                table = {},
                num_records = len(table),
                config = config,
                table_name = subtable_name,
                )
            table.set_cache(self.module_name, subtable)

        name_in_module = self.output_map[name]
        self.resolver.resolve_column(name_in_module, subtable)

        futures = subtable.get_futures(name_in_module)
        futures = wrap_future_exceptions(futures, table.table_name, name)

        table.set_futures(column, futures)




####

        
# Base class for resolvers that calculate descriptors

class Calculator(Resolver):
    input_names = []  # All of the names must exist for a calculation to occur
    output_names = []

    # Dispatch function for the normal case where the inputs
    # must not be exceptions. Pass only the exception-free
    # results to the actual calculation function, and use an
    # API which is a bit easier to add incrementally.
    
    def resolve_column(self, name, table):
        column_name = name

        # Get all of the inputs and see which ones have realizable values
        input_futures = [table.get_futures(input_name)
                                for input_name in self.input_names]
        input_values = []
        valid_record_indices = []
        n = len(table)

        # All of the results, as futures, will be placed here.
        result_columns = [[None]*n for name in self.output_names]

        for record_index, futures in enumerate(zip(*input_futures)):
            try:
                # Get the actual values
                values = []
                for future_index, future in enumerate(futures):
                    values.append(future.result())
            except Exception, err:
                # Couldn't get a value.
                # Mark all of the calculation children as "could not compute"
                failure_name = self.input_names[future_index]
                for output_name, result_column in zip(self.output_names, result_columns):
                    result_column[record_index] = _new_future_exception(
                        ResolverError(future.exception(), table.table_name, output_name))
            else:
                # Got a value. Add that to the list of terms to compute
                # and keep track of the mapping from the to-compute list
                # to the result columns
                valid_record_indices.append(record_index)
                input_values.append(values)

        if valid_record_indices:
            # There's something to calculate. Prepare the OutputDescriptors,
            # which is a support object for the calculation function, then
            # call the function.
            valid_record_indices.reverse()  # reverse so I can pop
            output_descriptors = OutputDescriptors(
                column_name, table.table_name, valid_record_indices, self.output_names, result_columns)

            self.calculate(table, name, input_values, output_descriptors)

        if valid_record_indices:
            # Looks like the calculation code didn't go to completion.
            # This shouldn't happen.
            raise AssertionError("incomplete processing: %d for %r"
                                 % (len(result_columns), descriptor))
        
        # Copy the results to the main table.
        for output_name, futures in zip(self.output_names, result_columns):
            table.set_futures(output_name, futures)

# This is a support object to make it easier to calculate descriptors.

class OutputDescriptors(object):
    def __init__(self, column, table_name, valid_record_indices, output_names, result_columns):
        self._column = column
        self._table_name = table_name
        self._valid_record_indices = valid_record_indices
        self._output_names = output_names
        self._result_columns = result_columns
        self._n = len(output_names)

    def add_result(self, result):
        """Specify the result for the next record.

        Can only be used when the calculator returns a single value.
        """
        if self._n != 1:
            raise TypeError("Can only use add_result() when there is a single output value")
        next_index = self._valid_record_indices.pop()
        self._result_columns[0][next_index] = _new_future(result)

    def add_results(self, results):
        """Specify a list/tuple of result values for the next record.

        The results must be in the same order as the output_names for the calculator.
        """
        if len(results) != self._n:
            raise TypeError("Expected %d descriptors but got %d"
                            % (self._n, len(results)))
        next_index = self._valid_record_indices.pop()
        for column, result in zip(self._result_columns, results):
            column[next_index] = _new_future(result)
            
    def add_futures(self, futures):
        """Specify a list/tuple of futures for the next record.

        The futures must be in the same order as the output_names for the calculator.
        """
        if len(futures) != self._n:
            raise TypeError("Expected %d future descriptors but got %d"
                            % (self._n, len(futures)))
        next_index = self._valid_record_indices.pop()
        for column, future in zip(self._result_columns, futures):
            column[next_index] = _new_future(future)

    def add_exception(self, exception):
        """Specify an exception for the next record."""
        next_index = self._valid_record_indices.pop()
        for output_name, column in zip(self._output_names, self._result_columns):
            column[next_index] = _new_future_exception(
                ResolverError(exception, self._table_name, output_name))
                #CalculationError(exception, output_name, self._table_name))

    def add_column_results(self, column_results):
        """Specify a list/tuple of values for each column

        Use this to specify all of the values for the remaining records.
        Each element of column_results must contain a list of size N,
        where N is the number of records in the request. The first list
        contains all of the values for the descriptor output_names[0],
        etc. The value is wrapped in a future and placed into the table.
        """
        if len(column_results) != len(self._output_names):
            raise ValueError("Expecting %d items in column_results, got %d"
                            % (len(self._output_names), len(column_results)))
        n = len(self._valid_record_indices)
        for i, column in enumerate(column_results):
            if len(column) != n:
                raise TypeError("Expecting %d values for column %r, got %d"
                                % (n, self._output_names[i], len(column)))

        # un-reverse the indices of where to place the results
        indices = self._valid_record_indices[::-1]
        # Set up the mapping from valid record index to results index
        pairs = list(enumerate(indices))
        # Put all of the futures in the right location
        for column, result_column in zip(column_results, self._result_columns):
            for i, index in pairs:
                result_column[index] = _new_future(column[i])
        
    def add_column_futures(self, column_futures):
        """Specify a list/tuple of futures for each column

        Use this to specify all of the futures for the remaining records.
        Each element of column_futures must contain a list of size N,
        where N is the number of records in the request. The first list
        contains all of the values for the descriptor output_names[0], etc.
        """
        if len(column_futures) != len(self._output_names):
            raise ValueError("Expecting %d items in column_futures, got %d"
                            % (len(self._output_names), len(column_futures)))
        n = len(self._valid_record_indices)
        for i, column in enumerate(columns):
            if len(column) != n:
                raise TypeError("Expecting %d values for column %r, got %d"
                                % (n, self._output_names[i], len(column)))

        # un-reverse the indices of where to place the results
        indices = self._valid_record_indices[::-1]
        # Set up the mapping from valid record index to results index
        pairs = list(enumerate(indices))
        # Put all of the futures in the right location
        for column, result_column in zip(column_futures, self._result_columns):
            for i, index in pairs:
                result_column[index] = column[i]

# A calculator which calls a function that computes a single descriptor name
                
class CalculateName(Calculator):
    def __init__(self, input_names, output_name, f, docstring, include_table=False,
                 resolver_name=None):
        self.input_names = input_names
        self.output_names = [output_name]
        self.f = f
        self.docstring = docstring
        self.include_table = include_table
        self.resolver_name = resolver_name

    def __repr__(self):
        return "CalculateName(%r, %r, %r, %r, %r, %r)" % (
            self.input_names, self.output_names[0], self.f, self.docstring,
            self.include_table, self.resolver_name)
    
    def calculate(self, table, name, input_values, output_descriptors):
        f = self.f
        include_table = self.include_table
        for values in input_values:
            if include_table:
                values += (table,)
            try:
                output_descriptors.add_results([f(*values)])
            except Exception, err:
                if DEBUG:
                    import traceback, sys
                    sys.stderr.write("Unable to calculate %r using %r with arguments %r\n"
                                     % (table.get_qualified_name(name), f, tuple(values)))
                    traceback.print_exc()
                output_descriptors.add_exception(err)


# A calculator which calls a function that computes multiple descriptor names
    
class CalculateNames(Calculator):
    def __init__(self, input_names, output_names, f, docstring, include_table=False,
                 resolver_name=None):
        self.input_names = input_names
        self.output_names = output_names
        self.f = f
        self.docstring = docstring
        self.include_table = include_table
        self.resolver_name = resolver_name

    def __repr__(self):
        return "CalculateNames(%r, %r, %r, %r, %r, %r)" % (
            self.input_names, self.output_names, self.f, self.docstring,
            self.include_table, self.resolver_name)
    
    def calculate(self, table, name, input_values, output_descriptors):
        f = self.f
        include_table = self.include_table
        for values in input_values:
            if include_table:
                values += (table,)
            try:
                output_descriptors.add_results(f(*values))
            except Exception, err:
                if DEBUG:
                    import traceback, sys
                    sys.stderr.write("Unable to calculate %r using %r with arguments %r\n"
                                     % (table.get_qualified_name(name), f, tuple(values)))
                    traceback.print_exc()
                output_descriptors.add_exception(err)
        

def calculate(input_names=None, output_names=None,
              include_table=False, docstring=None):
    """Decorator to indicate that a function can be used as a calculator resolver.

    The function will be called once for each input record, for those
    records where the input future does not contain an exception.
              
    `input_names` describes the columns used to get the input values.
    If it's list of names then the function will be called with the
    correponding descriptor values, in the same order as input_names. As
    a shorthand, a single string is the same as a one-element list
    containing that string. If input_names is None then the function
    argument names will be used as the list of input names, though if
    `include_table` is True then the last parameter name will be ignored.


    `output_names` is either a string or a list. If it's a string then
    the function's return value will be stored under that name, and the
    function will be wrapped inside of a CalculateName instance. If
    `output_names` a list of names then the function's return value must
    contain the same number of values, and each value is stored under its
    associated name.

    `include_table` specifies if the property table should be included
    as the last argument in the request. This is useful to get
    configuration information from the table.
    
    `docstring` is a description of the calculation function. It isn't
    yet used. In the future I want some way to query a resolver and
    report the individual steps.
    """
    if not (include_table is True or include_table is False):
        raise ValueError("include_table must be True or False: %r" % (include_table,))
    
    if input_names is not None:
        if isinstance(input_names, basestring):
            input_names = [input_names]
        elif isinstance(input_names, tuple):
            pass
        elif isinstance(input_names, list):
            input_names = tuple(input_names)
        else:
            raise ValueError("input_names must be None, a string, or a list/tuple, not %r"
                             % (input_names,))

    if output_names is not None:
        if isinstance(output_names, basestring):
            pass
        elif isinstance(output_names, tuple):
            pass
        elif isinstance(output_names, list):
            output_names = tuple(output_names)
        else:
            raise ValueError("output_names must be None, a string, or a list/tuple, not %r"
                             % (output_names,))
    
    def decorate_propbox_function(f):
        import inspect
        if input_names is None:
            argspec = inspect.getargspec(f)
            input_names_ = argspec.args
            if include_table:
                if not input_names_:
                    raise ValueError("requested 'include_table' but no parameter defined for it")
                input_names_ = input_names_[:-1]
        else:
            input_names_ = input_names

        if output_names is None:
            f_name = f.__name__
            if f_name[:5] == "calc_":
                if len(f_name) == 5:
                    raise ValueError(
                        "decorator does not specify output_names and the function name %r is too short"
                        % (f_name,))
                output_names_ = f_name[5:]
            else:
                raise ValueError(
                    "decorator does not specify output_names and the function name %r does not start with 'calc_'"
                    % (f_name,))
        else:
            output_names_ = output_names

        if docstring is None:
            docstring_ = inspect.getdoc(f)
        else:
            docstring_ = docstring

        resolver_name = getattr(f, "__name__", None)
        
        if isinstance(output_names_, basestring):
            calc = CalculateName(input_names_, output_names_, f, docstring, include_table, resolver_name)
        else:
            calc = CalculateNames(input_names_, output_names_, f, docstring, include_table, resolver_name)
        f.propbox_resolver = calc
        return f
    
    return decorate_propbox_function


def collect_resolvers(d=None):
    if d is None:
        import inspect
        d = inspect.stack()[1][0].f_globals
    resolvers = []
    for name, value in d.iteritems():
        # OEChem places some odd things in the local namespace.
        # This will be fixed in OpenEye's mid-2015 release.
        # Until then, ignore odd things.
        if not isinstance(name, basestring):
            continue

        if isinstance(value, Resolver):
            resolver = value
        else:
            resolver = getattr(value, "propbox_resolver", None)
            if resolver is None:
                continue
        resolvers.append(resolver)
    return Propbox(resolvers)

##########

_delimiters = {
    "tab": "\t",
    "space": " ",
    "whitespace": " ",
    }                                                          

class PropertyTable(object):
    def __init__(self, resolver, table, num_records, config=None, table_name=None):
        self.resolver = resolver
        self.table = table
        self._num_records = num_records
        if config is None:
            config = {}
        self.config = config
        self.table_name = table_name

        self._cache = {}

    def __len__(self):
        return self._num_records

    def get_qualified_name(self, name):
        if self.table_name is None:
            return name
        return "%s.%s" % (self.table_name, name)

    def get_futures(self, name):
        if name in self.table:
            return self.table[name]
            
        self.resolver.resolve_column(name, self)
        try:
            return self.table[name]
        except KeyError:
            raise AssertionError("Resolver %r did not set values for column %r"
                                 % (self.resolver, name))

    def set_futures(self, column, futures):
        if column in self.table:
            raise ValueError("Column name %r already set" % (column,))
        if len(futures) != self._num_records:
            raise ValueError("Column requires a list of size %d but got %d elements"
                             % (self._num_records, len(futures)))
        self.table[column] = futures
            
    def get_values(self, name, error_value=None):
        futures = self.get_futures(name)
        values = []
        for future in futures:
            if future.exception() is not None:
                values.append(error_value)
            else:
                values.append(future.result())
        return values

    def set_values(self, name, values):
        if name in self.table:
            raise ValueError("Column name %r already set" % (name,))
        if len(values) != self._num_records:
            raise ValueError("Column requires a list of size %d but got %d elements"
                             % (self._num_records, len(values)))
        self.table[name] = [_new_future(value) for value in values]
        

    def get_future_records(self, names):
        table_columns = []
        for name in names:
            if name not in self.table:
                self.resolver.resolve_column(name, self)
            try:
                table_columns.append(self.table[name])
            except KeyError:
                raise AssertionError(
                    "Resolver %r did not set descriptor for column %r"
                    % (self.resolver, name))

        records = []
        for values in zip(*table_columns):
            record = dict(zip(names, values))
            records.append(record)
        return records


    def get_records(self, names, exception_value=None):
        new_records = []
        for record in self.get_future_records(names):
            new_record = {}
            for k, v in record.iteritems():
                if v.exception() is not None:
                    new_record[k] = exception_value
                else:
                    new_record[k] = v.result()
            new_records.append(new_record)
        return new_records

    def add_future_records(self, future_records):
        if len(future_records) != self._num_records:
            raise AssertionError("Table contains %d records but trying to set %d"
                                 % (self._num_records), len(future_records))
        if not future_records:
            return
        
        # Only keys from the first record will be used
        # If other records contain additional keys, they will be ignored
        keys = set(future_records[0])
        for key in keys:
            if key in self.table:
                raise AssertionError("Table already contains the descriptor name %r"
                                     % (key,))
        
        d = {key: [] for key in keys} 
        try:
            for i, future_record in enumerate(future_records):
                for key in keys:
                    d[key].append(future_record[key])
        except KeyError:
            raise AssertionError("Record %d is missing descriptor name %r. It contains: %r"
                                 % (i, key, sorted(future_record)))
                
        for key in keys:
            self.table[key] = d[key]

    def has_computed(self, name):
        return name in self.table

    def get_cache_value(self, key):
        return self._cache.get(key, None)

    def set_cache_value(self, key, value):
        self._cache[key] = value

    def save(self, destination, names, headers=None,
             missing_values=None, formatters=None, include_header=True, dialect="tab"):
        if headers is None:
            headers = names
        elif len(names) != len(headers):
            raise ValueError("'headers' contains %d elements, but there are %d column names"
                             % (len(headers), len(names)))


        # Get the formatter function for each
        if formatters is None:
            formatters = [str] * len(names)
        elif len(formatters) != len(names):
            raise ValueError("'formatters' contains %d elements, but there are %d column names"
                             % (len(formatters), len(names)))
                
        # Get the string to use when something is missing
        if missing_values is None:
            missing_values = ["*"] * len(names)
        elif len(missing_values) != len(names):
            raise ValueError("'missing_values' contains %d elements, but there are %d column names"
                             % (len(missing_values), len(names)))
        else:
            for missing_value in missing_values:
                if not isinstance(missing_value, basestring):
                    raise ValueError("'missing_values' contains the non-string %r" % (missing_value,))

        # Get the columns to display
        missing_obj = object()
        columns = [self.get_values(name, missing_obj) for name in names]

        # Convert them to a string, or use the appropriate missing value
        new_columns = []
        for (name, column, formatter, missing_value) in zip(name, columns, formatters, missing_values):
            values = []
            try:
                for row_id, value in enumerate(column):
                    if value is missing_obj:
                        value = missing_value
                    else:
                        value = formatter(value)
                    if not isinstance(value, str):
                        raise ValueError("formatter did not return a string")
                    values.append(value)
            except Exception, err:
                raise ValueError("Cannot format column %r row %d (%r) as a %r: %s"
                                 % (name, row_id, value, formatter, err))
            new_columns.append(values)


        # When given a filename, I need to open and close it.
        if isinstance(destination, basestring):
            fileobj = open(destination, "w")
            close = fileobj.close
        else:
            fileobj = destination
            close = None

        try:
            # See if it's one of the easy delimiters
            delimiter = _delimiters.get(dialect, None)
            if delimiter is not None:
                # It is! I can write it myself.
                if include_header:
                    header = delimiter.join(headers) + "\n"
                    fileobj.write(header)
                for values in zip(*new_columns):
                    line = delimiter.join(values) + "\n"
                    fileobj.write(line)
            
            else:
                # Use the csv module
                import csv
                writer = csv.writer(fileobj, dialect)

                if include_header:
                    writer.writerow(headers)
                for values in zip(*new_columns):
                    writer.writerow(values)
        finally:
            if close is not None:
                close()
        
####
        
        
def _records_to_columns(records):
    # The records come in as a dictionary per compound.
    # I store the data as columns, with all values for a given descriptor.
    # This also verifies that all record dictionaries have the 
    # same descriptor terms.
    if not records:
        return {}
    columns = None
    n = 0
    for record_i, record in enumerate(records):
        if columns is None:
            columns = dict((descriptor, [value]) for descriptor, value in record.items())
            n = len(columns)
        else:
            if len(record) != n:
                raise ValueError("Record %d has %d elements (%r), expecting %d (%r)"
                                 % (record_i, len(record), sorted(record),
                                    n, sorted(columns)))
            try:
                for descriptor, value in record.items():
                    columns[descriptor].append(value)
            except KeyError:
                raise ValueError("Record %d contains the unexpected descriptor %r"
                                 % (record_i, descriptor))
    return columns

def _import_resolver(resolver):
    if not isinstance(resolver, basestring):
        return resolver
    module_resolver = __import__(resolver, globals(), {}, ["resolver"], 0)
    return module_resolver

def _add_id_column(future_columns, num_records):
    if "id" not in future_columns:
        future_columns["id"] = [_new_future("ID%d" % i) for i in range(1, num_records+1)]
    

def make_table_from_records(resolver, initial_records, config=None):
    resolver = _import_resolver(resolver)
    num_records = len(initial_records)
    columns = _records_to_columns(initial_records)
    future_columns = {}
    for descriptor, values in columns.iteritems():
        future_columns[descriptor] = [_new_future(value) for value in values]
    _add_id_column(future_columns, num_records)

    return PropertyTable(resolver, future_columns, num_records,
                         config=config)

def make_table_from_columns(resolver, initial_columns, config=None):
    resolver = _import_resolver(resolver)
    values = initial_columns.values()
    if not values:
        # Short-circuit the empty table.
        return PropertyTable(resolver, {}, 0, config=config)
    
    # Double-check that all of the list lengths are the same
    lengths = set(map(len, values))
    if len(lengths) == 1:
        future_columns = {}
        for descriptor, values in initial_columns.iteritems():
            if hasattr(values[0], "add_done_callback"):
                future_columns[descriptor] = values # Is this too hackish?
            else:
                future_columns[descriptor] = [_new_future(value) for value in values]
        num_records = list(lengths)[0]
        _add_id_column(future_columns, num_records)
        
        return PropertyTable(resolver, future_columns, num_records, config=config)

    # Report the problem column
        
    n = None
    reference = None
    # Sort so it's in a consistent order
    for descriptor, values in sorted(initial_columns.items()):
        if n is None:
            n = len(values)
            reference = descriptor
        else:
            if n != len(values):
                raise ValueError("Column length mismatch: %r has length %d while %r has length %d"
                                 % (reference, n,
                                    descriptor, len(values)))
    raise AssertionError("record list had different lengths but now they are the same")
