from datetime import timedelta

BROKER_BACKEND = "redis"
BROKER_URL = "redis://localhost:6379/0"
CELERY_RESULT_BACKEND = "redis://localhost:6379/0"

CELERY_TRACK_STARTED = True
CELERY_TASK_SERIALIZER = "json"
CELERY_RESULT_SERIALIZER = "json"

CELERYCAM_EXPIRE_SUCCESS = timedelta(weeks=520)
CELERYCAM_EXPIRE_ERROR = timedelta(weeks=520)
CELERYCAM_EXPIRE_PENDING = timedelta(weeks=520)

BROKER_TRANSPORT_OPTIONS = {'visibility_timeout': 172800} # 2 days

#CELERY_ALWAYS_EAGER = True
#CELERY_EAGER_PROPAGATES_EXCEPTIONS = True

#CELERY_IMPORTS = ("talent.tasks", )
