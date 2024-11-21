- ### Q: Since the nucleotides are case-insensitive, do we have to remember the case of the nucleotides? Or can we always show uppercase letters, but accept both uppercase and lowercase letters as input?
  - #### A: Input can be anything, but output can be what we want

- ### Q: Is the qkmer type's purpose to be stored in the DB or is it only used for querying?
  - #### A: Only querying, but still need to create a type

- ### Q: The cast for the DNA datatype seems to be implicit, do we need to add the cast explicitly ? 
  - ### A: We should add an explicit cast

- ### Q: Error management, do we need to return null or doing it like we did id okay ? (We did it the same way as in complex)
  - ### A: We should terminate on an error
