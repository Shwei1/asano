class InvalidPolygonError(Exception):
    """
    Помилка для неправильно введеного багатокутника
    """
    def __init__(self, message):
        self.message = message
        super().__init__(message)
